//! Threading helpers (std only), and the CPU name for reports.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock};

pub fn num_threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

/// Threads worth starting for `work` simple operations: starting a thread
/// costs about as much as 20,000 of them.
pub fn threads_for(work: usize, max: usize) -> usize {
    (work / 20_000).clamp(1, max.max(1))
}

/// Runs `f(offset, chunk)` over pieces of `data` on `threads` threads.
/// Pieces are handed out as threads finish, because cores differ in speed
/// (Apple silicon has performance and efficiency cores).
pub fn for_chunks<T, F>(data: &mut [T], threads: usize, f: F)
where
    T: Send,
    F: Fn(usize, &mut [T]) + Sync,
{
    let n = data.len();
    if threads <= 1 || n < 2 {
        f(0, data);
        return;
    }
    let size = n.div_ceil(threads * 8);
    let queue = Mutex::new(data.chunks_mut(size).enumerate());
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let next = queue.lock().unwrap().next();
                match next {
                    Some((k, chunk)) => f(k * size, chunk),
                    None => break,
                }
            });
        }
    });
}

/// Runs `f(item, &mut state)` for items 0..n on `threads` threads, each with
/// its own state from `init`, handing out items as threads finish.
/// Returns the states.
pub fn for_items<S, I, F>(n: usize, threads: usize, init: I, f: F) -> Vec<S>
where
    S: Send,
    I: Fn() -> S + Sync,
    F: Fn(usize, &mut S) + Sync,
{
    let threads = threads.clamp(1, n.max(1));
    let next = AtomicUsize::new(0);
    let work = || {
        let mut state = init();
        loop {
            let i = next.fetch_add(1, Ordering::Relaxed);
            if i >= n {
                return state;
            }
            f(i, &mut state);
        }
    };
    if threads == 1 {
        return vec![work()];
    }
    std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads).map(|_| s.spawn(work)).collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    })
}

pub fn cpu_name() -> &'static str {
    static NAME: OnceLock<String> = OnceLock::new();
    NAME.get_or_init(|| {
        std::process::Command::new("sysctl")
            .args(["-n", "machdep.cpu.brand_string"])
            .output()
            .ok()
            .and_then(|o| String::from_utf8(o.stdout).ok())
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "this CPU".to_string())
    })
}
