//! Fixed-point distance matrix for the in-memory engine.

use rayon::prelude::*;

use super::DistanceCalculator;
use crate::io::phylip::{PhylipMatrix, SCALE};

/// The strict upper triangle of a symmetric distance matrix, stored as
/// fixed-point `i32` values in units of `1e-8`, rounded to a multiple of 100.
///
/// Row `i` holds the distances to `i+1..k`, packed contiguously. The
/// in-memory engine mutates entries in place as nodes are merged.
///
/// # Examples
///
/// ```
/// use ninja::distance::{DistanceCalculator, DistanceMatrix};
/// use ninja::io::fasta::parse_fasta;
///
/// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n>c\nAGGT\n", None).unwrap();
/// let calc = DistanceCalculator::new(&aln, None).unwrap();
/// let m = DistanceMatrix::from_calculator(&calc);
/// assert_eq!(m.get(0, 1), m.get(1, 0));
/// assert_eq!(m.len(), 3);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DistanceMatrix {
    k: usize,
    data: Vec<i32>,
}

impl DistanceMatrix {
    /// Fixed-point scale of the stored values.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// assert_eq!(DistanceMatrix::SCALE, 100_000_000);
    /// assert_eq!(DistanceMatrix::quantize(1.0), DistanceMatrix::SCALE as i32);
    /// ```
    pub const SCALE: i64 = SCALE;

    /// Convert a real-valued distance to the stored fixed-point form.
    ///
    /// This is `100 * trunc((d * 1e8 + 50) / 100)`: round half up to the
    /// nearest `1e-6`, which is the precision the Phylip writer prints, so a
    /// tree built from an alignment and one built from the written matrix
    /// see identical distances.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// assert_eq!(DistanceMatrix::quantize(0.25), 25_000_000);
    /// assert_eq!(DistanceMatrix::quantize(0.123_456_785), 12_345_700);
    /// ```
    #[inline]
    pub fn quantize(d: f64) -> i32 {
        100 * (((SCALE as f64 * d) + 50.0) / 100.0) as i32
    }

    /// An all-zero matrix for `k` taxa.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// let m = DistanceMatrix::zeros(3);
    /// assert_eq!(m.len(), 3);
    /// assert_eq!(m.get(0, 1), 0);
    /// ```
    pub fn zeros(k: usize) -> Self {
        DistanceMatrix { k, data: vec![0; Self::tri_len(k)] }
    }

    fn tri_len(k: usize) -> usize {
        k * k.saturating_sub(1) / 2
    }

    /// Compute every pairwise distance from an alignment, in parallel.
    ///
    /// Rows are distributed across the current rayon thread pool.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::{DistanceCalculator, DistanceMatrix};
    /// use ninja::io::fasta::parse_fasta;
    ///
    /// let aln = parse_fasta(b">a\nACGT\n>b\nACGA\n>c\nAGGT\n", None).unwrap();
    /// let calc = DistanceCalculator::new(&aln, None).unwrap();
    /// let m = DistanceMatrix::from_calculator(&calc);
    /// assert_eq!(m.as_slice().len(), 3);
    /// ```
    pub fn from_calculator(calc: &DistanceCalculator) -> Self {
        let k = calc.len();
        let mut m = Self::zeros(k);
        // Split `data` into per-row slices of decreasing length so rows can
        // be filled independently.
        let mut rows: Vec<&mut [i32]> = Vec::with_capacity(k);
        let mut rest: &mut [i32] = &mut m.data;
        for i in 0..k {
            let (row, tail) = rest.split_at_mut(k - i - 1);
            rows.push(row);
            rest = tail;
        }
        rows.into_par_iter().enumerate().for_each(|(i, row)| {
            for (off, cell) in row.iter_mut().enumerate() {
                let j = i + 1 + off;
                *cell = Self::quantize(calc.calc(i, j));
            }
        });
        m
    }

    /// Build from a parsed Phylip matrix (already in `1e-8` units).
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    /// use ninja::io::phylip::read_phylip_from;
    ///
    /// let lower = "3\na\nb 0.5\nc 0.25 0.75\n";
    /// let p = read_phylip_from(lower.as_bytes()).unwrap();
    /// let m = DistanceMatrix::from_phylip(&p);
    /// assert!((m.get_f64(0, 2) - 0.25).abs() < 1e-9);
    /// ```
    pub fn from_phylip(p: &PhylipMatrix) -> Self {
        let k = p.len();
        let mut m = Self::zeros(k);
        for i in 0..k {
            for j in (i + 1)..k {
                let v = p.lower[j][i];
                // Match the quantisation applied to computed distances.
                let v = 100 * ((v + 50) / 100);
                m.data[Self::index(k, i, j)] = v.clamp(i32::MIN as i64, i32::MAX as i64) as i32;
            }
        }
        m
    }

    /// Number of taxa.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// assert_eq!(DistanceMatrix::zeros(5).len(), 5);
    /// ```
    #[inline]
    pub fn len(&self) -> usize {
        self.k
    }

    /// True when there are no taxa.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// assert!(DistanceMatrix::zeros(0).is_empty());
    /// assert!(!DistanceMatrix::zeros(2).is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.k == 0
    }

    #[inline]
    fn index(k: usize, i: usize, j: usize) -> usize {
        debug_assert!(i < j && j < k);
        i * (2 * k - i - 1) / 2 + (j - i - 1)
    }

    /// Hint the cache that `(i, j)` will be read soon.
    #[inline]
    #[allow(unsafe_code)]
    pub(crate) fn prefetch(&self, i: usize, j: usize) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        let idx = Self::index(self.k, a, b);
        #[cfg(target_arch = "x86_64")]
        {
            let p = self.data.as_ptr().wrapping_add(idx) as *const i8;
            // SAFETY: prefetch has no architectural effect and takes any
            // address; the pointer is derived from a live slice.
            unsafe { std::arch::x86_64::_mm_prefetch(p, std::arch::x86_64::_MM_HINT_T0) };
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = idx;
        }
    }

    /// Distance between `i` and `j` (`i != j`), in fixed-point units.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// let mut m = DistanceMatrix::zeros(3);
    /// m.set(0, 1, 500);
    /// assert_eq!(m.get(0, 1), 500);
    /// assert_eq!(m.get(1, 0), 500);
    /// ```
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> i32 {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.data[Self::index(self.k, a, b)]
    }

    /// Overwrite the distance between `i` and `j` (`i != j`).
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// let mut m = DistanceMatrix::zeros(3);
    /// m.set(1, 2, 42);
    /// assert_eq!(m.get(2, 1), 42);
    /// ```
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: i32) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        let k = self.k;
        self.data[Self::index(k, a, b)] = v;
    }

    /// Distance between `i` and `j` as a real number.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// let mut m = DistanceMatrix::zeros(3);
    /// m.set(0, 1, DistanceMatrix::quantize(0.25));
    /// assert!((m.get_f64(0, 1) - 0.25).abs() < 1e-9);
    /// ```
    pub fn get_f64(&self, i: usize, j: usize) -> f64 {
        self.get(i, j) as f64 / SCALE as f64
    }

    /// The packed upper triangle.
    ///
    /// # Examples
    ///
    /// ```
    /// use ninja::distance::DistanceMatrix;
    ///
    /// let mut m = DistanceMatrix::zeros(3);
    /// m.set(0, 1, 7);
    /// assert_eq!(m.as_slice(), &[7, 0, 0]);
    /// ```
    pub fn as_slice(&self) -> &[i32] {
        &self.data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantize_rounds_to_hundreds() {
        assert_eq!(DistanceMatrix::quantize(0.705596), 70_559_600);
        assert_eq!(DistanceMatrix::quantize(0.7055964), 70_559_600);
        assert_eq!(DistanceMatrix::quantize(0.7055965), 70_559_700);
        assert_eq!(DistanceMatrix::quantize(0.0), 0);
        assert_eq!(DistanceMatrix::quantize(3.0), 300_000_000);
    }

    #[test]
    fn indexing_round_trips() {
        let k = 7;
        let mut m = DistanceMatrix::zeros(k);
        let mut v = 1;
        for i in 0..k {
            for j in (i + 1)..k {
                m.set(i, j, v);
                v += 1;
            }
        }
        assert_eq!(m.as_slice().len(), 21);
        let mut v = 1;
        for i in 0..k {
            for j in (i + 1)..k {
                assert_eq!(m.get(i, j), v);
                assert_eq!(m.get(j, i), v);
                v += 1;
            }
        }
    }
}
