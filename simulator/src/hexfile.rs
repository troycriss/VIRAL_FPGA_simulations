//! The hex files the Verilog design reads: INIT_FILE and WEIGHT_FILE.

/// One hex value per line, as $readmemh reads them (blank lines and
/// // comments are skipped).
pub fn read(path: &str) -> Result<Vec<u64>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut values = Vec::new();
    for (k, line) in text.lines().enumerate() {
        let line = line.split("//").next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let v = u64::from_str_radix(line, 16)
            .map_err(|_| format!("{path}:{}: not a hex number: {line}", k + 1))?;
        values.push(v);
    }
    Ok(values)
}

/// INIT_FILE: {x, y, z} per line, each a 16-bit Q1.15 value.
pub fn read_spins(path: &str) -> Result<Vec<[i16; 3]>, String> {
    Ok(read(path)?
        .into_iter()
        .map(|v| [(v >> 32) as u16 as i16, (v >> 16) as u16 as i16, v as u16 as i16])
        .collect())
}

/// WEIGHT_FILE: unsigned, 15 fractional bits (8000 = 1.0).
pub fn read_weights(path: &str) -> Result<Vec<u16>, String> {
    read(path)?
        .into_iter()
        .map(|v| u16::try_from(v).map_err(|_| format!("{path}: weight {v:x} is wider than 16 bits")))
        .collect()
}

/// A spin as written by the testbench's %h: 12 lowercase hex digits.
pub fn spin_hex(s: [i16; 3]) -> String {
    format!("{:04x}{:04x}{:04x}", s[0] as u16, s[1] as u16, s[2] as u16)
}
