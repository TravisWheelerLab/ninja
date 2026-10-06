//! FASTA alignment reader.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use crate::alphabet::Alphabet;
use crate::error::{Error, Result};

/// A multiple sequence alignment: one row per sequence, all the same length.
///
/// # Examples
///
/// ```
/// use ninja::io::fasta::parse_fasta;
///
/// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n", None).unwrap();
/// assert_eq!(aln.names, vec!["a", "b"]);
/// assert_eq!(aln.seqs, vec![b"ACGT".to_vec(), b"ACGA".to_vec()]);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// Sequence identifiers: the text after `>` up to the first whitespace.
    pub names: Vec<String>,
    /// Upper-case residues, one `Vec<u8>` per sequence. Gaps are `-`.
    pub seqs: Vec<Vec<u8>>,
    /// Alphabet, either supplied or detected from the residues.
    pub alphabet: Alphabet,
}

impl Alignment {
    /// Number of sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n", None).unwrap();
    /// assert_eq!(aln.len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// True when the alignment holds no sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::alphabet::Alphabet;
    /// use ninja::io::fasta::{parse_fasta, Alignment};
    ///
    /// let empty = Alignment { names: Vec::new(), seqs: Vec::new(), alphabet: Alphabet::Dna };
    /// assert!(empty.is_empty());
    ///
    /// let aln = parse_fasta(b">a\nACGT\n", None).unwrap();
    /// assert!(!aln.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Number of columns (after all-gap columns were removed).
    ///
    /// # Examples
    ///
    /// Columns 0 and 1 are gaps in every sequence, so they are dropped and
    /// `width` reports 2, not the 4 columns in the input.
    ///
    /// ```
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\n--AC\n>b\n--GT\n", None).unwrap();
    /// assert_eq!(aln.width(), 2);
    /// assert_eq!(aln.seqs[0], b"AC");
    /// ```
    pub fn width(&self) -> usize {
        self.seqs.first().map_or(0, |s| s.len())
    }
}

/// Read a FASTA alignment from a file.
///
/// See [`read_fasta_from`] for the parsing rules.
///
/// # Examples
///
/// ```
/// use std::io::Write;
///
/// use ninja::io::fasta::read_fasta;
///
/// let mut file = tempfile::NamedTempFile::new().unwrap();
/// write!(file, ">a\nACGT\n>b\nACGA\n").unwrap();
/// let aln = read_fasta(file.path(), None).unwrap();
/// assert_eq!(aln.names, vec!["a", "b"]);
/// ```
pub fn read_fasta(path: impl AsRef<Path>, alphabet: Option<Alphabet>) -> Result<Alignment> {
    let path = path.as_ref();
    let f = File::open(path).map_err(|e| Error::io(path, e))?;
    let mut bytes = Vec::new();
    BufReader::new(f).read_to_end(&mut bytes).map_err(|e| Error::io(path, e))?;
    parse_fasta(&bytes, alphabet)
}

/// Read a FASTA alignment from any reader (for example standard input).
///
/// # Examples
///
/// ```
/// use std::io::Cursor;
///
/// use ninja::io::fasta::read_fasta_from;
///
/// let aln = read_fasta_from(Cursor::new(b">a\nACGT\n>b\nACGA\n"), None).unwrap();
/// assert_eq!(aln.names, vec!["a", "b"]);
/// ```
pub fn read_fasta_from(mut reader: impl Read, alphabet: Option<Alphabet>) -> Result<Alignment> {
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes)?;
    parse_fasta(&bytes, alphabet)
}

/// Parse FASTA text.
///
/// Rules, following the reference implementation:
///
/// * A record starts at `>`; the name is everything up to the first space,
///   tab, or newline. The rest of the header line is ignored.
/// * Every non-whitespace byte on the following lines is a residue. `.` is
///   treated as a gap and stored as `-`.
/// * Residues are converted to upper case.
/// * Columns that are gaps in every sequence are removed.
/// * If no alphabet is given, the alignment is DNA when every non-gap residue
///   is one of `A C G T U`, otherwise protein. (The reference checked case
///   before upper-casing, so lower-case DNA was mistaken for protein; this
///   reader checks after upper-casing.)
/// * For DNA, `U` is converted to `T`.
///
/// All sequences must have the same length, and there must be at least one.
///
/// # Examples
///
/// The header's description is dropped, `.` and `-` both mean a gap,
/// lower-case residues are upper-cased, the detected DNA alphabet turns the
/// trailing `U` into `T`, and the column that ends up a gap in both
/// sequences is removed.
///
/// ```
/// use ninja::alphabet::Alphabet;
/// use ninja::io::fasta::parse_fasta;
///
/// let aln = parse_fasta(b">a desc\nAC-G\nT\n>b\nac.gu\n", None).unwrap();
/// assert_eq!(aln.names, vec!["a", "b"]);
/// assert_eq!(aln.alphabet, Alphabet::Dna);
/// assert_eq!(aln.seqs, vec![b"ACGT".to_vec(), b"ACGT".to_vec()]);
/// ```
pub fn parse_fasta(bytes: &[u8], alphabet: Option<Alphabet>) -> Result<Alignment> {
    let mut names: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();

    let mut i = 0;
    let n = bytes.len();
    // Skip anything before the first record.
    while i < n && bytes[i] != b'>' {
        if !bytes[i].is_ascii_whitespace() {
            return Err(Error::format(
                "expected FASTA input starting with '>' (only FASTA alignments are accepted)",
            ));
        }
        i += 1;
    }
    while i < n {
        debug_assert_eq!(bytes[i], b'>');
        i += 1;
        let name_start = i;
        while i < n && !matches!(bytes[i], b' ' | b'\t' | b'\n' | b'\r') {
            i += 1;
        }
        let name = String::from_utf8_lossy(&bytes[name_start..i]).into_owned();
        if name.is_empty() {
            return Err(Error::format(format!("empty sequence name in record {}", names.len() + 1)));
        }
        // Skip the remainder of the header line.
        while i < n && bytes[i] != b'\n' {
            i += 1;
        }
        let mut seq = Vec::new();
        while i < n && bytes[i] != b'>' {
            let c = bytes[i];
            if !c.is_ascii_whitespace() {
                seq.push(if c == b'.' { b'-' } else { c.to_ascii_uppercase() });
            }
            i += 1;
        }
        names.push(name);
        seqs.push(seq);
    }

    if names.is_empty() {
        return Err(Error::invalid("no sequences found in input"));
    }
    let width = seqs[0].len();
    for (k, s) in seqs.iter().enumerate() {
        if s.len() != width {
            return Err(Error::invalid(format!(
                "sequence '{}' has length {} but '{}' has length {}; input must be an alignment",
                names[k],
                s.len(),
                names[0],
                width
            )));
        }
    }

    // Drop all-gap columns.
    let keep: Vec<bool> = (0..width).map(|c| seqs.iter().any(|s| s[c] != b'-')).collect();
    if keep.iter().any(|&k| !k) {
        for s in seqs.iter_mut() {
            let mut w = 0;
            for c in 0..width {
                if keep[c] {
                    s[w] = s[c];
                    w += 1;
                }
            }
            s.truncate(w);
        }
    }

    let alphabet = alphabet.unwrap_or_else(|| detect_alphabet(&seqs));
    if alphabet == Alphabet::Dna {
        for s in seqs.iter_mut() {
            for c in s.iter_mut() {
                if *c == b'U' {
                    *c = b'T';
                }
            }
        }
    }

    Ok(Alignment { names, seqs, alphabet })
}

impl Alignment {
    /// Group sequences with identical residues (after all-gap columns were
    /// removed), each group in input order and listed by its first member.
    ///
    /// # Examples
    ///
    /// `a`, `c`, and `d` share a sequence, so they form one group; `b` forms
    /// its own group.
    ///
    /// ```
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n>c\nACGT\n>d\nACGT\n", None).unwrap();
    /// assert_eq!(aln.duplicate_groups(), vec![vec![0, 2, 3], vec![1]]);
    /// ```
    pub fn duplicate_groups(&self) -> Vec<Vec<usize>> {
        use std::collections::HashMap;
        let mut first: HashMap<&[u8], usize> = HashMap::new();
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for (i, s) in self.seqs.iter().enumerate() {
            match first.get(s.as_slice()) {
                Some(&g) => groups[g].push(i),
                None => {
                    first.insert(s.as_slice(), groups.len());
                    groups.push(vec![i]);
                }
            }
        }
        groups
    }

    /// The alignment restricted to the first member of each group.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n>c\nACGT\n>d\nACGT\n", None).unwrap();
    /// let reps = aln.representatives(&aln.duplicate_groups());
    /// assert_eq!(reps.names, vec!["a", "b"]);
    /// assert_eq!(reps.seqs, vec![b"ACGT".to_vec(), b"ACGA".to_vec()]);
    /// ```
    pub fn representatives(&self, groups: &[Vec<usize>]) -> Alignment {
        Alignment {
            names: groups.iter().map(|g| self.names[g[0]].clone()).collect(),
            seqs: groups.iter().map(|g| self.seqs[g[0]].clone()).collect(),
            alphabet: self.alphabet,
        }
    }
}

/// DNA when every non-gap residue is in `ACGTU`, otherwise protein.
///
/// # Examples
///
/// ```
/// use ninja::alphabet::Alphabet;
/// use ninja::io::fasta::detect_alphabet;
///
/// assert_eq!(detect_alphabet(&[b"ACGT".to_vec(), b"AC-U".to_vec()]), Alphabet::Dna);
/// assert_eq!(detect_alphabet(&[b"ACDE".to_vec()]), Alphabet::Amino);
/// ```
pub fn detect_alphabet(seqs: &[Vec<u8>]) -> Alphabet {
    let dna = seqs.iter().all(|s| s.iter().all(|&c| matches!(c, b'-' | b'A' | b'C' | b'G' | b'T' | b'U')));
    if dna {
        Alphabet::Dna
    } else {
        Alphabet::Amino
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_names_and_strips_gap_columns() {
        let txt = b">s1 desc\nAC-G\nT\n>s2\nac.gu\n";
        let a = parse_fasta(txt, None).unwrap();
        assert_eq!(a.names, vec!["s1", "s2"]);
        assert_eq!(a.alphabet, Alphabet::Dna);
        assert_eq!(a.seqs, vec![b"ACGT".to_vec(), b"ACGT".to_vec()]);
    }

    #[test]
    fn detects_protein() {
        let a = parse_fasta(b">a\nACDE\n>b\nACDF\n", None).unwrap();
        assert_eq!(a.alphabet, Alphabet::Amino);
    }

    #[test]
    fn rejects_ragged() {
        assert!(parse_fasta(b">a\nACDE\n>b\nACD\n", None).is_err());
    }

    #[test]
    fn rejects_empty() {
        assert!(parse_fasta(b"", None).is_err());
        assert!(parse_fasta(b"ACGT\n", None).is_err());
    }

    #[test]
    fn duplicate_groups() {
        let a = parse_fasta(b">a\nACGT\n>b\nACGA\n>c\nACGT\n>d\nACGT\n", None).unwrap();
        assert_eq!(a.duplicate_groups(), vec![vec![0, 2, 3], vec![1]]);
        let r = a.representatives(&a.duplicate_groups());
        assert_eq!(r.names, vec!["a", "b"]);
    }
}
