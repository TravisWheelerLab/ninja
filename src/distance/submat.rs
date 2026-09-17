//! Amino-acid substitution matrices and the dissimilarities derived from them.
//!
//! A log-odds score `s(a, b)` in units of `bits` bits encodes the odds ratio
//! `q(a, b) / (p(a) p(b))`, so its reciprocal `2^(-s·bits)` is a
//! dissimilarity. Dividing that by the geometric mean of the two residues'
//! average dissimilarity to a random residue, under a fixed background
//! composition, puts the expected dissimilarity of a random pair near 1.
//! FastTree derived its BLOSUM45 table this way; the derivation here
//! reproduces that table to 1e-14 before rounding to `f32`.

use std::path::Path;

use crate::error::{Error, Result};

/// Residues in table order.
pub const RESIDUES: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

/// Background composition per mille, in [`RESIDUES`] order: the average
/// protein composition that FastTree's table implies (Doolittle's table,
/// with glutamine and glutamate in the order it was entered there).
const BACKGROUND_PERMILLE: [u32; 20] =
    [74, 42, 44, 59, 33, 58, 37, 74, 29, 38, 76, 72, 18, 40, 50, 81, 62, 13, 33, 68];

/// A 20x20 integer scoring matrix and its score unit.
#[derive(Debug, Clone, PartialEq)]
pub struct SubstitutionMatrix {
    name: String,
    scores: [[i32; 20]; 20],
    bits: f64,
    scale_declared: bool,
    builtin_blosum45: bool,
}

impl SubstitutionMatrix {
    /// BLOSUM62 in half-bit units, the default for protein distances.
    pub fn blosum62() -> Self {
        Self::parse(include_str!("matrices/BLOSUM62.txt"), "BLOSUM62").expect("built-in BLOSUM62 parses")
    }

    /// BLOSUM45 in third-bit units. Its dissimilarities are FastTree's
    /// published table rather than a fresh derivation, so output matches
    /// earlier releases bit for bit.
    pub fn blosum45() -> Self {
        let mut m =
            Self::parse(include_str!("matrices/BLOSUM45.txt"), "BLOSUM45").expect("built-in BLOSUM45 parses");
        m.builtin_blosum45 = true;
        m
    }

    /// A built-in matrix by name, case-insensitively: `BLOSUM62` or `BLOSUM45`.
    pub fn by_name(name: &str) -> Option<Self> {
        match name.to_ascii_uppercase().as_str() {
            "BLOSUM62" => Some(Self::blosum62()),
            "BLOSUM45" => Some(Self::blosum45()),
            _ => None,
        }
    }

    /// Read a matrix file in the NCBI format used by BLAST and EMBOSS.
    pub fn from_file(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|e| Error::io(path, e))?;
        let name =
            path.file_name().map_or_else(|| path.display().to_string(), |f| f.to_string_lossy().into_owned());
        Self::parse(&text, &name)
    }

    /// Parse the NCBI matrix format: `#` comment lines, a header row of
    /// residue letters, then one row per residue with an integer score per
    /// column. Columns and rows for `B`, `Z`, `X` and `*` are ignored. The
    /// score unit is taken from a comment such as `scale of ln(2)/2` or
    /// `in 1/2 Bit Units`; without one it is estimated from the scores and
    /// [`scale_declared`](Self::scale_declared) reports `false`.
    pub fn parse(text: &str, name: &str) -> Result<Self> {
        let mut bits: Option<f64> = None;
        let mut columns: Option<Vec<Option<usize>>> = None;
        let mut scores = [[i32::MIN; 20]; 20];
        let mut has_row = [false; 20];
        for (k, raw) in text.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(comment) = line.strip_prefix('#') {
                if bits.is_none() {
                    bits = scale_from_comment(comment);
                }
                continue;
            }
            let bad = |msg: String| Error::format(format!("{}: line {}: {}", name, k + 1, msg));
            let fields: Vec<&str> = line.split_whitespace().collect();
            match &columns {
                None => {
                    let mut map = Vec::with_capacity(fields.len());
                    for f in &fields {
                        if f.len() != 1 {
                            return Err(bad(format!(
                                "expected residue letters in the header row, found '{}'",
                                f
                            )));
                        }
                        map.push(residue_index(f.as_bytes()[0]));
                    }
                    columns = Some(map);
                }
                Some(map) => {
                    let (row, values) = fields.split_first().expect("non-empty line has a field");
                    if row.len() != 1 {
                        return Err(bad(format!(
                            "expected a residue letter at the start of the row, found '{}'",
                            row
                        )));
                    }
                    if values.len() != map.len() {
                        return Err(bad(format!("expected {} scores, found {}", map.len(), values.len())));
                    }
                    let Some(a) = residue_index(row.as_bytes()[0]) else { continue };
                    has_row[a] = true;
                    for (col, v) in map.iter().zip(values) {
                        let Some(b) = *col else { continue };
                        scores[a][b] =
                            v.parse().map_err(|_| bad(format!("'{}' is not an integer score", v)))?;
                    }
                }
            }
        }
        let has_col = |b: usize| columns.as_ref().is_some_and(|map| map.contains(&Some(b)));
        let missing: Vec<String> = (0..20)
            .filter(|&a| !has_row[a] || !has_col(a))
            .map(|a| (RESIDUES[a] as char).to_string())
            .collect();
        if !missing.is_empty() {
            return Err(Error::format(format!(
                "{}: no scores for residue(s) {} (a matrix needs a row and a column for each of {})",
                name,
                missing.join(", "),
                std::str::from_utf8(RESIDUES).unwrap()
            )));
        }
        let (bits, scale_declared) = match bits {
            Some(b) => (b, true),
            None => (estimate_bits(&scores), false),
        };
        Ok(SubstitutionMatrix {
            name: name.to_string(),
            scores,
            bits,
            scale_declared,
            builtin_blosum45: false,
        })
    }

    /// The name given when the matrix was created (a built-in name or a file name).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Bits per score unit (1/2 for BLOSUM62, 1/3 for BLOSUM45).
    pub fn bits(&self) -> f64 {
        self.bits
    }

    /// Whether the score unit came from the file rather than an estimate.
    pub fn scale_declared(&self) -> bool {
        self.scale_declared
    }

    /// Score for residues `a` and `b`, indices into [`RESIDUES`].
    pub fn score(&self, a: usize, b: usize) -> i32 {
        self.scores[a][b]
    }

    /// The dissimilarity table used for distances: FastTree's published
    /// table for the built-in BLOSUM45, the derivation otherwise.
    pub fn dissimilarities(&self) -> [[f32; 20]; 20] {
        if self.builtin_blosum45 {
            super::bl45::BL45
        } else {
            self.derived_dissimilarities()
        }
    }

    /// The derivation described in the module documentation, for any matrix.
    pub fn derived_dissimilarities(&self) -> [[f32; 20]; 20] {
        let freq: Vec<f64> = BACKGROUND_PERMILLE.iter().map(|&p| p as f64 / 1001.0).collect();
        let mut m = [[0f64; 20]; 20];
        for a in 0..20 {
            for b in 0..20 {
                if a != b {
                    m[a][b] = 2f64.powf(-(self.scores[a][b] as f64) * self.bits);
                }
            }
        }
        let r: Vec<f64> = (0..20).map(|a| (0..20).map(|b| freq[b] * m[a][b]).sum()).collect();
        let mut d = [[0f32; 20]; 20];
        for a in 0..20 {
            for b in 0..20 {
                if a != b {
                    d[a][b] = (m[a][b] / (r[a] * r[b]).sqrt()) as f32;
                }
            }
        }
        d
    }
}

fn residue_index(letter: u8) -> Option<usize> {
    RESIDUES.iter().position(|&r| r == letter.to_ascii_uppercase())
}

/// Bits per score unit from a header comment: `... scale of ln(2)/3.0`
/// gives 1/3, `... in 1/2 Bit Units` gives 1/2, `... in Bit Units` gives 1.
fn scale_from_comment(comment: &str) -> Option<f64> {
    let lower = comment.to_ascii_lowercase();
    let number = |s: &str| -> Option<f64> {
        let n: String = s.trim_start().chars().take_while(|c| c.is_ascii_digit() || *c == '.').collect();
        n.trim_end_matches('.').parse::<f64>().ok().filter(|x| *x > 0.0)
    };
    if let Some(i) = lower.find("ln(2)/") {
        return number(&lower[i + 6..]).map(|x| 1.0 / x);
    }
    if let Some(i) = lower.find("bit unit") {
        let before = lower[..i].trim_end();
        if let Some(j) = before.rfind("1/") {
            return number(&before[j + 2..]).map(|x| 1.0 / x);
        }
        if before.ends_with(" in") {
            return Some(1.0);
        }
    }
    None
}

/// Solve `sum_ab f_a f_b 2^(lambda s_ab) = 1` for the nonzero root by
/// bisection, using the background composition above. An estimate only:
/// the matrix's own background would give the exact unit.
fn estimate_bits(scores: &[[i32; 20]; 20]) -> f64 {
    let freq: Vec<f64> = BACKGROUND_PERMILLE.iter().map(|&p| p as f64 / 1001.0).collect();
    let g = |lambda: f64| -> f64 {
        let mut total = 0.0;
        for a in 0..20 {
            for b in 0..20 {
                total += freq[a] * freq[b] * 2f64.powf(lambda * scores[a][b] as f64);
            }
        }
        total - 1.0
    };
    let mut lo = 1e-3;
    if g(lo) >= 0.0 {
        return 0.5;
    }
    let mut hi = 0.1;
    while g(hi) < 0.0 {
        hi *= 2.0;
        if hi > 64.0 {
            return 0.5;
        }
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if g(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::bl45::BL45;

    #[test]
    fn blosum45_derivation_reproduces_fasttree_table() {
        let m = SubstitutionMatrix::blosum45();
        assert_eq!(m.bits(), 1.0 / 3.0);
        let d = m.derived_dissimilarities();
        let mut worst = 0f32;
        for a in 0..20 {
            for b in 0..20 {
                worst = worst.max((d[a][b] - BL45[a][b]).abs());
            }
        }
        assert!(worst < 2e-6, "max difference {worst}");
        assert_eq!(m.dissimilarities(), BL45);
    }

    #[test]
    fn blosum62_parses_with_its_scale() {
        let m = SubstitutionMatrix::blosum62();
        assert!(m.scale_declared());
        assert_eq!(m.bits(), 0.5);
        let w = residue_index(b'W').unwrap();
        assert_eq!(m.score(w, w), 11);
        assert_eq!(m.score(0, 1), -1);
        assert_eq!(m.score(1, 0), -1);
        let d = m.dissimilarities();
        assert_eq!(d[0][0], 0.0);
        assert!(d[0][1] > 0.5 && d[0][1] < 3.0, "{}", d[0][1]);
        assert!(SubstitutionMatrix::by_name("blosum62").is_some());
        assert!(SubstitutionMatrix::by_name("PAM250").is_none());
    }

    #[test]
    fn scale_is_estimated_when_not_declared() {
        let text: String = include_str!("matrices/BLOSUM45.txt")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n");
        let m = SubstitutionMatrix::parse(&text, "test").unwrap();
        assert!(!m.scale_declared());
        assert!((m.bits() - 1.0 / 3.0).abs() < 0.02, "{}", m.bits());
    }

    #[test]
    fn scale_comments_are_understood() {
        assert_eq!(
            scale_from_comment("  Entries for the BLOSUM62 matrix at a scale of ln(2)/2.0."),
            Some(0.5)
        );
        assert_eq!(scale_from_comment("  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units"), Some(1.0 / 3.0));
        assert_eq!(scale_from_comment("  scores in Bit Units"), Some(1.0));
        assert_eq!(scale_from_comment("  Blocks Database = blocks.dat"), None);
    }

    #[test]
    fn missing_residue_is_an_error() {
        let text = "   A  R\nA  4 -1\nR -1  5\n";
        let err = SubstitutionMatrix::parse(text, "small").unwrap_err().to_string();
        assert!(err.contains("no scores for residue(s) N, D"), "{err}");
        let err = SubstitutionMatrix::parse("   A  R\nA  4\n", "short").unwrap_err().to_string();
        assert!(err.contains("line 2") && err.contains("expected 2 scores"), "{err}");
    }
}
