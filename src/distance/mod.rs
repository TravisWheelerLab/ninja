//! Pairwise evolutionary distances from an alignment.
//!
//! The calculator packs each sequence once so that a pair distance is a tight
//! loop over machine words:
//!
//! * DNA is stored two bits per site (`A=00, G=01, C=10, T=11`) plus a
//!   one-bit-per-site validity mask, 32 sites to a `u64`. XOR of two packed
//!   sequences yields `01` exactly for transitions (A<->G, C<->T) and a set
//!   high bit for transversions, so transitions, transversions, and the
//!   number of comparable sites are three population counts. This replaces
//!   the hand-written SSE shuffle kernel of the C port with code that
//!   vectorises on any target.
//! * Protein is stored one byte per site as an index into a 21x21 BLOSUM45
//!   dissimilarity table whose row and column 20 (any non-standard residue)
//!   are zero, so the inner loop is a branch-free gather and accumulate.
//!
//! The arithmetic (single-precision accumulation, the order of promotions in
//! the correction formulas) follows the reference implementation so that
//! distances agree with it to the printed precision.

mod bl45;
mod dna;
mod gaps;
mod matrix;
mod protein;
pub mod submat;

pub use matrix::DistanceMatrix;
pub use submat::SubstitutionMatrix;

use crate::alphabet::{Alphabet, Correction};
use crate::error::{Error, Result};
use crate::io::fasta::Alignment;

/// Computes corrected distances between any two sequences of an alignment.
///
/// # Examples
///
/// ```
/// use ninja::distance::DistanceCalculator;
/// use ninja::io::fasta::parse_fasta;
///
/// let aln = parse_fasta(b">a\nACGTACGTAC\n>b\nACGTACGTAG\n>c\nAGCTACGTAC\n", None).unwrap();
/// let calc = DistanceCalculator::new(&aln, None).unwrap();
/// let d = calc.calc(0, 1);
/// assert!((d - 0.108_466_1).abs() < 1e-6);
/// ```
#[derive(Debug, Clone)]
pub struct DistanceCalculator {
    alphabet: Alphabet,
    correction: Correction,
    packed: Packed,
    n: usize,
}

#[derive(Debug, Clone)]
enum Packed {
    Dna(dna::PackedDna),
    Amino(protein::PackedProtein),
}

impl DistanceCalculator {
    /// Pack an alignment for distance computation.
    ///
    /// `correction` defaults to Kimura two-parameter for DNA and scoredist
    /// for protein. A correction that does not apply to the alphabet is an
    /// error.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::alphabet::Correction;
    /// use ninja::distance::DistanceCalculator;
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGTACGTAC\n>b\nACGTACGTAG\n>c\nAGCTACGTAC\n", None).unwrap();
    /// let calc = DistanceCalculator::new(&aln, Some(Correction::JukesCantor)).unwrap();
    /// assert!((calc.calc(0, 1) - 0.107_325_6).abs() < 1e-6);
    /// ```
    pub fn new(aln: &Alignment, correction: Option<Correction>) -> Result<Self> {
        Self::with_matrix(aln, correction, &SubstitutionMatrix::blosum62())
    }

    /// Like [`new`](Self::new) with a chosen substitution matrix for
    /// protein distances (ignored for DNA).
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::{DistanceCalculator, SubstitutionMatrix};
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nMKVLATIR\n>b\nMKVLSTVR\n>c\nMRVLATIK\n", None).unwrap();
    /// let calc = DistanceCalculator::with_matrix(&aln, None, &SubstitutionMatrix::blosum45()).unwrap();
    /// assert!((calc.calc(0, 1) - 0.183_226_1).abs() < 1e-6);
    /// ```
    pub fn with_matrix(
        aln: &Alignment,
        correction: Option<Correction>,
        matrix: &SubstitutionMatrix,
    ) -> Result<Self> {
        let alphabet = aln.alphabet;
        let correction = correction.unwrap_or_else(|| Correction::default_for(alphabet));
        if !correction.applies_to(alphabet) {
            return Err(Error::options(format!(
                "correction '{}' cannot be used with the {} alphabet",
                correction, alphabet
            )));
        }
        let packed = match alphabet {
            Alphabet::Dna => Packed::Dna(dna::PackedDna::new(&aln.seqs)),
            Alphabet::Amino => {
                Packed::Amino(protein::PackedProtein::new(&aln.seqs, &matrix.dissimilarities()))
            }
        };
        Ok(DistanceCalculator { alphabet, correction, packed, n: aln.len() })
    }

    /// Number of sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceCalculator;
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n>c\nAGGT\n", None).unwrap();
    /// let calc = DistanceCalculator::new(&aln, None).unwrap();
    /// assert_eq!(calc.len(), 3);
    /// ```
    pub fn len(&self) -> usize {
        self.n
    }

    /// True when there are no sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::alphabet::Alphabet;
    /// use ninja::distance::DistanceCalculator;
    /// use ninja::io::fasta::Alignment;
    ///
    /// let aln = Alignment { names: vec![], seqs: vec![], alphabet: Alphabet::Dna };
    /// let calc = DistanceCalculator::new(&aln, None).unwrap();
    /// assert!(calc.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The alphabet the calculator was built for.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::alphabet::Alphabet;
    /// use ninja::distance::DistanceCalculator;
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n", None).unwrap();
    /// let calc = DistanceCalculator::new(&aln, None).unwrap();
    /// assert_eq!(calc.alphabet(), Alphabet::Dna);
    /// ```
    pub fn alphabet(&self) -> Alphabet {
        self.alphabet
    }

    /// The correction in use.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::alphabet::Correction;
    /// use ninja::distance::DistanceCalculator;
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n", None).unwrap();
    /// let calc = DistanceCalculator::new(&aln, Some(Correction::None)).unwrap();
    /// assert_eq!(calc.correction(), Correction::None);
    /// ```
    pub fn correction(&self) -> Correction {
        self.correction
    }

    /// Corrected distance between sequences `a` and `b`.
    ///
    /// Pairs with no comparable sites get the correction's maximum distance
    /// (1 with no correction, 3 otherwise), as do pairs whose correction
    /// formula is undefined (saturated divergence).
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceCalculator;
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGTACGTAC\n>b\nACGTACGTAG\n>c\nAGCTACGTAC\n", None).unwrap();
    /// let calc = DistanceCalculator::new(&aln, None).unwrap();
    /// assert!((calc.calc(0, 2) - 0.239_278_2).abs() < 1e-6);
    /// ```
    #[inline]
    pub fn calc(&self, a: usize, b: usize) -> f64 {
        let maxscore = self.correction.max_distance();
        match (&self.packed, self.correction) {
            (Packed::Dna(p), Correction::OneGap) => {
                let (transitions, transversions, sites) = p.count(a, b);
                let openings = p.openings(a, b);
                gaps::onegap_distance(transitions + transversions, sites, openings, maxscore) as f64
            }
            (Packed::Amino(p), Correction::OneGap) => {
                let (mismatches, sites, openings) = p.count_onegap(a, b);
                gaps::onegap_distance(mismatches, sites, openings, maxscore) as f64
            }
            (Packed::Dna(p), corr) => {
                let (transitions, transversions, sites) = p.count(a, b);
                dna::correct(transitions, transversions, sites, corr) as f64
            }
            (Packed::Amino(p), corr) => {
                let (sum, sites) = p.score(a, b);
                protein::correct(sum, sites, corr) as f64
            }
        }
    }
}
