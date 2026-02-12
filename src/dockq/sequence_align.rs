//! Needleman-Wunsch global sequence alignment for chain matching.
//!
//! Provides sequence alignment to establish residue correspondence
//! between model and native chains in DockQ calculations.

/// Parameters for Needleman-Wunsch alignment scoring.
#[derive(Debug, Clone)]
pub struct AlignmentParams {
    /// Score for matching residues.
    pub match_score: i32,
    /// Score for mismatching residues.
    pub mismatch_score: i32,
    /// Penalty for opening a gap.
    pub gap_open: i32,
    /// Penalty for extending a gap.
    pub gap_extend: i32,
}

impl Default for AlignmentParams {
    fn default() -> Self {
        Self {
            match_score: 5,
            mismatch_score: 0,
            gap_open: -4,
            gap_extend: -1,
        }
    }
}

/// Result of a sequence alignment.
#[derive(Debug, Clone)]
pub struct SequenceAlignment {
    /// Mapping from seq1 positions to seq2 positions (None = gap).
    pub mapping: Vec<Option<usize>>,
    /// Alignment score.
    pub score: i32,
    /// Fraction of aligned positions that are identical.
    pub identity: f64,
    /// Number of aligned (non-gap) positions.
    pub num_aligned: usize,
}

/// Perform Needleman-Wunsch global alignment with affine gap penalties.
///
/// Returns an alignment mapping positions in `seq1` to positions in `seq2`.
pub fn align_sequences(
    seq1: &[String],
    seq2: &[String],
    params: &AlignmentParams,
) -> SequenceAlignment {
    let n = seq1.len();
    let m = seq2.len();

    if n == 0 || m == 0 {
        return SequenceAlignment {
            mapping: vec![None; n],
            score: 0,
            identity: 0.0,
            num_aligned: 0,
        };
    }

    // Score matrix and traceback
    // Using simple linear gap model for clarity (gap_open + gap_extend per gap position)
    let rows = n + 1;
    let cols = m + 1;
    let mut score = vec![vec![0i32; cols]; rows];
    let mut traceback = vec![vec![0u8; cols]; rows]; // 0=diag, 1=up, 2=left

    // Initialize first row and column
    for i in 1..rows {
        score[i][0] = params.gap_open + params.gap_extend * (i as i32 - 1);
        traceback[i][0] = 1; // up (gap in seq2)
    }
    for j in 1..cols {
        score[0][j] = params.gap_open + params.gap_extend * (j as i32 - 1);
        traceback[0][j] = 2; // left (gap in seq1)
    }

    // Fill matrix
    for i in 1..rows {
        for j in 1..cols {
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                params.match_score
            } else {
                params.mismatch_score
            };

            let diag = score[i - 1][j - 1] + match_mismatch;

            let up_gap_penalty = if traceback[i - 1][j] == 1 {
                params.gap_extend
            } else {
                params.gap_open
            };
            let up = score[i - 1][j] + up_gap_penalty;

            let left_gap_penalty = if traceback[i][j - 1] == 2 {
                params.gap_extend
            } else {
                params.gap_open
            };
            let left = score[i][j - 1] + left_gap_penalty;

            if diag >= up && diag >= left {
                score[i][j] = diag;
                traceback[i][j] = 0;
            } else if up >= left {
                score[i][j] = up;
                traceback[i][j] = 1;
            } else {
                score[i][j] = left;
                traceback[i][j] = 2;
            }
        }
    }

    // Traceback
    let mut mapping = vec![None; n];
    let mut i = n;
    let mut j = m;
    let mut num_identical = 0;
    let mut num_aligned = 0;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 && traceback[i][j] == 0 {
            // Diagonal: aligned pair
            mapping[i - 1] = Some(j - 1);
            num_aligned += 1;
            if seq1[i - 1] == seq2[j - 1] {
                num_identical += 1;
            }
            i -= 1;
            j -= 1;
        } else if i > 0 && traceback[i][j] == 1 {
            // Up: gap in seq2
            i -= 1;
        } else {
            // Left: gap in seq1
            j -= 1;
        }
    }

    let identity = if num_aligned > 0 {
        num_identical as f64 / num_aligned as f64
    } else {
        0.0
    };

    SequenceAlignment {
        mapping,
        score: score[n][m],
        identity,
        num_aligned,
    }
}

/// Compute sequence identity between two residue name sequences.
///
/// Uses Needleman-Wunsch alignment and returns the fraction of identical residues.
pub fn sequence_identity(seq1: &[String], seq2: &[String]) -> f64 {
    let alignment = align_sequences(seq1, seq2, &AlignmentParams::default());
    alignment.identity
}

#[cfg(test)]
mod tests {
    use super::*;

    fn to_strings(s: &[&str]) -> Vec<String> {
        s.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn test_identical_sequences() {
        let seq = to_strings(&["ALA", "GLY", "VAL", "LEU"]);
        let result = align_sequences(&seq, &seq, &AlignmentParams::default());

        assert_eq!(result.identity, 1.0);
        assert_eq!(result.num_aligned, 4);
        for (i, m) in result.mapping.iter().enumerate() {
            assert_eq!(*m, Some(i));
        }
    }

    #[test]
    fn test_completely_different() {
        let seq1 = to_strings(&["ALA", "ALA", "ALA"]);
        let seq2 = to_strings(&["GLY", "GLY", "GLY"]);
        let result = align_sequences(&seq1, &seq2, &AlignmentParams::default());

        assert_eq!(result.identity, 0.0);
    }

    #[test]
    fn test_one_mismatch() {
        let seq1 = to_strings(&["ALA", "GLY", "VAL"]);
        let seq2 = to_strings(&["ALA", "LEU", "VAL"]);
        let result = align_sequences(&seq1, &seq2, &AlignmentParams::default());

        // Should align position-by-position with one mismatch
        assert!(result.identity > 0.5);
        assert!(result.identity < 1.0);
    }

    #[test]
    fn test_different_lengths() {
        let seq1 = to_strings(&["ALA", "GLY", "VAL", "LEU", "ILE"]);
        let seq2 = to_strings(&["ALA", "GLY", "VAL"]);
        let result = align_sequences(&seq1, &seq2, &AlignmentParams::default());

        // Should still find the common subsequence
        assert!(result.num_aligned >= 3);
    }

    #[test]
    fn test_empty_sequence() {
        let seq1 = to_strings(&["ALA", "GLY"]);
        let seq2: Vec<String> = vec![];
        let result = align_sequences(&seq1, &seq2, &AlignmentParams::default());

        assert_eq!(result.identity, 0.0);
        assert_eq!(result.num_aligned, 0);
    }

    #[test]
    fn test_sequence_identity() {
        let seq1 = to_strings(&["ALA", "GLY", "VAL"]);
        let seq2 = to_strings(&["ALA", "GLY", "VAL"]);
        assert!((sequence_identity(&seq1, &seq2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_insertion_alignment() {
        // seq2 has an insertion relative to seq1
        let seq1 = to_strings(&["ALA", "GLY", "VAL"]);
        let seq2 = to_strings(&["ALA", "LEU", "GLY", "VAL"]);
        let result = align_sequences(&seq1, &seq2, &AlignmentParams::default());

        // ALA should map to ALA, GLY to GLY, VAL to VAL
        assert_eq!(result.mapping[0], Some(0)); // ALA -> ALA
        assert!(result.num_aligned >= 3);
    }
}
