use std::path::PathBuf;
use std::sync::atomic::AtomicBool;
use viral_sim::{bench, hexfile, hw, par, run, server, study};

const USAGE: &str = "\
viral-sim: CPU simulator for the VIRAL mean-field QAOA FPGA design

  viral-sim serve [--port 8080] [--web DIR]
      web page at http://127.0.0.1:8080/ (DIR: the repository's web/)

  viral-sim run [--KEY VALUE ...] [--json FILE]
      one run: prints a summary, and writes the web page's JSON with --json
      problem   partition (rank 1, runs on VIRAL_n_spin), sk, sk-gauss (dense),
                maxcut (random regular graph), chain (J = 1/r^alpha),
                hopfield (stored patterns, rank P)
      n         number of spins
      layers    number of layers
      dt        ramp size: gamma_k = s_k dt, beta_k = (1 - s_k) dt
      engine    hw (bit-exact VIRAL_n_spin, partition only) or float
      schedule  ramp, or lut (the hardware's 4-entry angle_lut)
      init      x (along +x, z tilted by up to +-tilt) or random
      degree    maxcut: edges per vertex (3)
      alpha     chain: J = 1/r^alpha (3)
      patterns, noise  hopfield: stored patterns (10), fraction of the cue flipped (0.2)
      also: seed, tilt, gamma-lo, points, threads

  viral-sim study --study scaling|sweep|squeezing [--KEY VALUE ...] [--json FILE]
      many runs, with results as curves:
      scaling    SK energy per spin against N, extrapolated to N -> infinity
                 (--sizes 250,500,...  --seeds 10, and the run options)
      sweep      final result over a grid of layers and dt
                 (--layers-list 25,50,...  --dt-list 0.1,0.2,...  --seeds 5)
      squeezing  one-axis twisting: many classical runs from random starts
                 (discrete truncated Wigner) against the exact quantum result
                 (--n 1000  --trajectories 2000  --times 60  --t-max)

  viral-sim hw --init FILE [--weights FILE] [--layers L] [--gamma-lo X]
               [--led-lsb 8] [--dump FILE]
      bit-exact VIRAL_n_spin run on the testbench's hex files, with the
      angle_lut schedule; --dump writes the final spins like the testbench

  viral-sim bench [--threads T]
      times layer updates on this CPU against the FPGA estimates
";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (cmd, rest): (&str, &[String]) = match args.split_first() {
        Some((c, r)) => (c.as_str(), r),
        None => ("help", &[]),
    };
    let result = match cmd {
        "serve" => serve(rest),
        "run" => run_cmd(rest),
        "hw" => hw_cmd(rest),
        "bench" => bench_cmd(rest),
        "study" => study_cmd(rest),
        _ => {
            print!("{USAGE}");
            Ok(())
        }
    };
    if let Err(e) = result {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

/// --key value ... -> [(key, value)], with - in keys turned into _.
fn pairs(args: &[String]) -> Result<Vec<(String, String)>, String> {
    let mut out = Vec::new();
    let mut it = args.iter();
    while let Some(k) = it.next() {
        let key = k.strip_prefix("--").ok_or_else(|| format!("expected --option, got {k:?}"))?;
        let v = it.next().ok_or_else(|| format!("--{key} needs a value"))?;
        out.push((key.replace('-', "_"), v.clone()));
    }
    Ok(out)
}

fn take(pairs: &mut Vec<(String, String)>, key: &str) -> Option<String> {
    let i = pairs.iter().position(|(k, _)| k == key)?;
    Some(pairs.remove(i).1)
}

fn parse<T: std::str::FromStr>(key: &str, v: Option<String>, default: T) -> Result<T, String> {
    match v {
        None => Ok(default),
        Some(v) => v.parse().map_err(|_| format!("--{key}: bad value {v:?}")),
    }
}

fn no_more(p: &[(String, String)]) -> Result<(), String> {
    match p.first() {
        Some((k, _)) => Err(format!("unknown option --{}", k.replace('_', "-"))),
        None => Ok(()),
    }
}

fn serve(args: &[String]) -> Result<(), String> {
    let mut p = pairs(args)?;
    let port = parse("port", take(&mut p, "port"), 8080u16)?;
    let web = take(&mut p, "web")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../web"));
    no_more(&p)?;
    let web = web.canonicalize().map_err(|e| format!("{}: {e}", web.display()))?;
    server::serve(port, web, par::num_threads())
}

fn run_cmd(args: &[String]) -> Result<(), String> {
    let mut p = pairs(args)?;
    let json = take(&mut p, "json");
    let params = run::Params::from_pairs(&p)?;
    let out = run::run(&params)?;
    println!("{}", out.summary);
    if let Some(path) = json {
        std::fs::write(&path, out.json.to_json()).map_err(|e| format!("{path}: {e}"))?;
        println!("wrote {path}");
    }
    Ok(())
}

fn hw_cmd(args: &[String]) -> Result<(), String> {
    let mut p = pairs(args)?;
    let init = take(&mut p, "init").ok_or("--init FILE is required")?;
    let weights = take(&mut p, "weights");
    let layers = parse("layers", take(&mut p, "layers"), 2usize)?;
    let gamma_lo: Option<u16> = take(&mut p, "gamma_lo")
        .map(|v| v.parse().map_err(|_| format!("--gamma-lo: bad value {v:?}")))
        .transpose()?;
    let led_lsb = parse("led-lsb", take(&mut p, "led_lsb"), 8u32)?;
    let dump = take(&mut p, "dump");
    no_more(&p)?;
    if layers > 4 {
        return Err("angle_lut has 4 entries: --layers at most 4".into());
    }

    let spins = hexfile::read_spins(&init)?;
    let n = spins.len();
    let a = match weights.as_deref() {
        Some(path) => {
            let w = hexfile::read_weights(path)?;
            if w.len() < n {
                return Err(format!("{path}: {} weights for {n} spins", w.len()));
            }
            w[..n].to_vec()
        }
        None => vec![0x8000; n],
    };
    let mut m = hw::Hw { spins, a };
    for k in 0..layers {
        let mut g = hw::LUT_GAMMA[k];
        if let Some(lo) = gamma_lo {
            g = (g & !0xFFFF) | lo as i32;
        }
        m.layer(g, hw::LUT_BETA[k], par::num_threads());
    }
    let clocks = hw::clocks(n as u64, layers as u64);
    let e = m.energy();
    let led = m.led(led_lsb);
    println!(
        "N_SPINS={n} N_LAYERS={layers}: {clocks} clocks, energy {e} ({:.6}), led {led:08b}",
        e as f64 / 32768.0
    );
    if let Some(path) = dump {
        let mut text: String = m.spins.iter().map(|&s| hexfile::spin_hex(s) + "\n").collect();
        text += &format!("energy {e}\nled {led:08b}\nclocks {clocks}\n");
        std::fs::write(&path, text).map_err(|e| format!("{path}: {e}"))?;
    }
    Ok(())
}

fn study_cmd(args: &[String]) -> Result<(), String> {
    let mut p = pairs(args)?;
    let json = take(&mut p, "json");
    let sp = study::StudyParams::from_pairs(&p)?;
    let report = |frac: f64, msg: &str, _: Option<viral_sim::json::Value>| {
        eprint!("\r{:5.1}%  {msg:<60}", 100.0 * frac);
    };
    let (mut v, text) = study::run(&sp, &report, &AtomicBool::new(false))?;
    if let viral_sim::json::Value::Obj(m) = &mut v {
        m.push(("summary".into(), viral_sim::json::Value::Str(text.clone())));
    }
    eprintln!();
    println!("{text}");
    if let Some(path) = json {
        std::fs::write(&path, v.to_json()).map_err(|e| format!("{path}: {e}"))?;
        println!("wrote {path}");
    }
    Ok(())
}

fn bench_cmd(args: &[String]) -> Result<(), String> {
    let mut p = pairs(args)?;
    let threads = parse("threads", take(&mut p, "threads"), par::num_threads())?;
    no_more(&p)?;
    bench::run(threads.max(1));
    Ok(())
}
