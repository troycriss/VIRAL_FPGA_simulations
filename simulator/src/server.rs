//! Local web server (std only): serves the page in web/ and runs
//! simulations for it at /api/run?problem=...&n=...
//!
//! Binds to 127.0.0.1 only, so nothing outside this computer can reach it.

use crate::json::Value;
use crate::run::{self, Params};
use crate::study::{self, StudyParams};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, Instant};

/// A study running in the background, polled by the page.
struct Job {
    cancel: AtomicBool,
    state: Mutex<JobState>,
}

#[derive(Default)]
struct JobState {
    progress: f64,
    message: String,
    /// latest partial result, serialized
    partial: Option<String>,
    /// final result or error
    result: Option<Result<String, String>>,
}

fn jobs() -> &'static Mutex<HashMap<u64, Arc<Job>>> {
    static JOBS: OnceLock<Mutex<HashMap<u64, Arc<Job>>>> = OnceLock::new();
    JOBS.get_or_init(|| Mutex::new(HashMap::new()))
}

static NEXT_JOB: AtomicU64 = AtomicU64::new(1);

/// /api/study/start, /status, /cancel
fn study_api(action: &str, pairs: &[(String, String)], threads: usize) -> Result<Value, String> {
    let id = || -> Result<u64, String> {
        pairs.iter().find(|(k, _)| k == "id").and_then(|(_, v)| v.parse().ok()).ok_or_else(|| "id: missing".to_string())
    };
    match action {
        "start" => {
            let mut pairs = pairs.to_vec();
            if !pairs.iter().any(|(k, _)| k == "threads") {
                pairs.push(("threads".into(), threads.to_string()));
            }
            let sp = StudyParams::from_pairs(&pairs)?;
            let job = Arc::new(Job { cancel: AtomicBool::new(false), state: Mutex::new(JobState::default()) });
            let id = NEXT_JOB.fetch_add(1, Ordering::Relaxed);
            {
                let mut all = jobs().lock().unwrap();
                // keep the 20 most recent jobs
                while all.len() >= 20 {
                    let oldest = *all.keys().min().unwrap();
                    all.remove(&oldest);
                }
                all.insert(id, job.clone());
            }
            std::thread::spawn(move || {
                let start = Instant::now();
                let report = |p: f64, msg: &str, partial: Option<Value>| {
                    let mut st = job.state.lock().unwrap();
                    st.progress = p;
                    st.message = msg.to_string();
                    if let Some(v) = partial {
                        st.partial = Some(v.to_json());
                    }
                };
                let r = study::run(&sp, &report, &job.cancel);
                println!("study {id}: {:.1} s", start.elapsed().as_secs_f64());
                let mut st = job.state.lock().unwrap();
                st.result = Some(r.map(|(v, text)| {
                    let mut v = v;
                    if let Value::Obj(m) = &mut v {
                        m.push(("summary".into(), Value::Str(text)));
                    }
                    v.to_json()
                }));
                if st.result.as_ref().is_some_and(|r| r.is_ok()) {
                    st.progress = 1.0;
                }
            });
            Ok(Value::obj().set("id", id))
        }
        "status" | "cancel" => {
            let id = id()?;
            let job = jobs().lock().unwrap().get(&id).cloned().ok_or("no such study (finished long ago?)")?;
            if action == "cancel" {
                job.cancel.store(true, Ordering::Relaxed);
            }
            let st = job.state.lock().unwrap();
            let mut v = Value::obj()
                .set("id", id)
                .set("progress", st.progress)
                .set("message", st.message.as_str())
                .set("done", st.result.is_some());
            match &st.result {
                Some(Ok(r)) => v = v.set("result", Value::Raw(r.clone())),
                Some(Err(e)) => v = v.set("error", e.as_str()),
                None => {
                    v = v.set("partial", st.partial.as_ref().map_or(Value::Null, |p| Value::Raw(p.clone())));
                }
            }
            Ok(v)
        }
        _ => Err(format!("unknown study action {action:?}")),
    }
}

pub fn serve(port: u16, web: PathBuf, threads: usize) -> Result<(), String> {
    let listener = TcpListener::bind(("127.0.0.1", port)).map_err(|e| format!("port {port}: {e}"))?;
    println!("Serving {} at http://127.0.0.1:{port}/  (Ctrl-C to stop)", web.display());
    for stream in listener.incoming() {
        let Ok(stream) = stream else { continue };
        let web = web.clone();
        std::thread::spawn(move || {
            if let Err(e) = handle(stream, &web, threads) {
                eprintln!("request failed: {e}");
            }
        });
    }
    Ok(())
}

fn handle(mut stream: TcpStream, web: &Path, threads: usize) -> std::io::Result<()> {
    stream.set_read_timeout(Some(Duration::from_secs(10)))?;
    let mut reader = BufReader::new(stream.try_clone()?);
    let mut request = String::new();
    reader.read_line(&mut request)?;
    loop {
        let mut header = String::new();
        if reader.read_line(&mut header)? == 0 || header.trim().is_empty() {
            break;
        }
    }
    let mut parts = request.split_whitespace();
    let (method, target) = (parts.next().unwrap_or(""), parts.next().unwrap_or("/"));
    if method != "GET" {
        return respond(&mut stream, 405, "text/plain", b"GET only");
    }
    let (path, query) = target.split_once('?').unwrap_or((target, ""));

    if path == "/api/run" {
        let start = Instant::now();
        let mut pairs = parse_query(query);
        if !pairs.iter().any(|(k, _)| k == "threads") {
            pairs.push(("threads".into(), threads.to_string()));
        }
        let result = Params::from_pairs(&pairs).and_then(|p| run::run(&p));
        println!("{target}: {:.2} s", start.elapsed().as_secs_f64());
        return match result {
            Ok(out) => respond(&mut stream, 200, "application/json", out.json.to_json().as_bytes()),
            Err(e) => {
                let body = Value::obj().set("error", e).to_json();
                respond(&mut stream, 400, "application/json", body.as_bytes())
            }
        };
    }

    if let Some(action) = path.strip_prefix("/api/study/") {
        let result = study_api(action, &parse_query(query), threads);
        if action == "start" {
            println!("{target}");
        }
        return match result {
            Ok(v) => respond(&mut stream, 200, "application/json", v.to_json().as_bytes()),
            Err(e) => {
                let body = Value::obj().set("error", e).to_json();
                respond(&mut stream, 400, "application/json", body.as_bytes())
            }
        };
    }

    // static files from web/, no paths outside it
    let rel = if path == "/" { "index.html" } else { path.trim_start_matches('/') };
    let safe = !rel.contains("..") && rel.chars().all(|c| c.is_ascii_alphanumeric() || "._-/".contains(c));
    match std::fs::read(web.join(rel)) {
        Ok(body) if safe => respond(&mut stream, 200, content_type(rel), &body),
        _ => respond(&mut stream, 404, "text/plain", b"not found"),
    }
}

fn content_type(path: &str) -> &'static str {
    match path.rsplit('.').next().unwrap_or("") {
        "html" => "text/html; charset=utf-8",
        "js" => "text/javascript; charset=utf-8",
        "css" => "text/css; charset=utf-8",
        "json" => "application/json",
        "svg" => "image/svg+xml",
        "png" => "image/png",
        _ => "application/octet-stream",
    }
}

fn respond(stream: &mut TcpStream, code: u16, ctype: &str, body: &[u8]) -> std::io::Result<()> {
    let reason = match code {
        200 => "OK",
        400 => "Bad Request",
        404 => "Not Found",
        _ => "Method Not Allowed",
    };
    let head = format!(
        "HTTP/1.1 {code} {reason}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\nCache-Control: no-store\r\nConnection: close\r\n\r\n",
        body.len()
    );
    stream.write_all(head.as_bytes())?;
    stream.write_all(body)?;
    stream.flush()
}

/// a=1&b=two%20words -> [(a, 1), (b, two words)]
pub fn parse_query(q: &str) -> Vec<(String, String)> {
    q.split('&')
        .filter(|s| !s.is_empty())
        .map(|kv| {
            let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
            (decode(k), decode(v))
        })
        .collect()
}

fn decode(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'+' => out.push(b' '),
            b'%' if i + 2 < b.len() => {
                match u8::from_str_radix(std::str::from_utf8(&b[i + 1..i + 3]).unwrap_or("zz"), 16) {
                    Ok(v) => {
                        out.push(v);
                        i += 2;
                    }
                    Err(_) => out.push(b'%'),
                }
            }
            c => out.push(c),
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}
