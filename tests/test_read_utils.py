"""
Tests for read_utils.py

Run with:  pytest tests/
"""
import numpy as np
import pytest

from read_utils import (
    annotate_positions,
    build_tally,
    classify_read_strand,
    is_chimeric,
    minimal_unique_prefix_length,
)


# ---------------------------------------------------------------------------
# build_tally
# ---------------------------------------------------------------------------

class TestBuildTally:
    def test_all_methylated(self):
        # lowercase XM chars = methylated
        tally = build_tally('zzz', [10, 11, 12], win_start=10, max_index=5)
        assert tally is not None
        np.testing.assert_array_equal(tally[:3], [1.0, 1.0, 1.0])
        assert np.all(np.isnan(tally[3:]))

    def test_all_unmethylated(self):
        tally = build_tally('ZZZ', [10, 11, 12], win_start=10, max_index=5)
        np.testing.assert_array_equal(tally[:3], [0.0, 0.0, 0.0])

    def test_mixed(self):
        tally = build_tally('ZzZ', [10, 11, 12], win_start=10, max_index=5)
        assert tally[0] == 0.0
        assert tally[1] == 1.0
        assert tally[2] == 0.0

    def test_non_alphabetic_chars_become_nan(self):
        # '.' is non-context (not alpha) → should be NaN
        tally = build_tally('Z.z', [10, 11, 12], win_start=10, max_index=5)
        assert tally[0] == 0.0
        assert np.isnan(tally[1])
        assert tally[2] == 1.0

    def test_positions_before_window_excluded(self):
        # positions [5,6,7] all have rel_pos < 0 relative to win_start=10
        result = build_tally('ZZZ', [5, 6, 7], win_start=10, max_index=5)
        assert result is None

    def test_positions_after_window_excluded(self):
        # only position 10 falls in [10, 15); 0 and 20 do not
        tally = build_tally('ZZZ', [0, 10, 20], win_start=10, max_index=5)
        assert tally[0] == 0.0         # position 10 → rel_pos 0
        assert np.all(np.isnan(tally[1:]))

    def test_no_valid_positions_returns_none(self):
        result = build_tally('ZZZ', [100, 200, 300], win_start=0, max_index=5)
        assert result is None

    def test_output_length_equals_max_index(self):
        tally = build_tally('Z', [0], win_start=0, max_index=10)
        assert tally.shape == (10,)

    def test_duplicate_positions_last_write_wins(self):
        # Two reads mapped to the same relative position: second overwrites first
        # 'Z' then 'z' at position 10 → methylated (1.0)
        tally = build_tally('Zz', [10, 10], win_start=10, max_index=3)
        assert tally[0] == 1.0


# ---------------------------------------------------------------------------
# classify_read_strand
# ---------------------------------------------------------------------------

class TestClassifyReadStrand:
    @pytest.mark.parametrize("motif,ys,expected", [
        ('+', 'OT',   'same'),
        ('+', 'CTOT', 'same'),
        ('+', 'CTOB', 'opp'),
        ('+', 'OB',   'opp'),
        ('-', 'OT',   'opp'),
        ('-', 'CTOT', 'opp'),
        ('-', 'CTOB', 'same'),
        ('-', 'OB',   'same'),
    ])
    def test_all_strand_combinations(self, motif, ys, expected):
        assert classify_read_strand(motif, ys) == expected

    def test_unknown_ys_tag_raises(self):
        with pytest.raises(ValueError):
            classify_read_strand('+', 'XX')

    def test_unknown_motif_strand_raises(self):
        with pytest.raises(ValueError):
            classify_read_strand('*', 'OT')

    def test_returns_string(self):
        result = classify_read_strand('+', 'OT')
        assert isinstance(result, str)


# ---------------------------------------------------------------------------
# annotate_positions
# ---------------------------------------------------------------------------

class TestAnnotatePositions:
    def _pairs(self, ref_seq, query_seq):
        """Build aligned_pairs for a gapless alignment."""
        return [(i, i, r) for i, r in enumerate(ref_seq)]

    def test_c_unconverted(self):
        arr = annotate_positions(self._pairs('C', 'C'), 'C', max_read_len=5)
        assert arr[0] == 'C'

    def test_c_converted_to_t(self):
        arr = annotate_positions(self._pairs('C', 'T'), 'T', max_read_len=5)
        assert arr[0] == 'T'

    def test_g_unconverted(self):
        arr = annotate_positions(self._pairs('G', 'G'), 'G', max_read_len=5)
        assert arr[0] == 'G'

    def test_g_converted_to_a(self):
        arr = annotate_positions(self._pairs('G', 'A'), 'A', max_read_len=5)
        assert arr[0] == 'A'

    def test_non_cg_match_is_dot(self):
        arr = annotate_positions(self._pairs('T', 'T'), 'T', max_read_len=5)
        assert arr[0] == '.'

    def test_non_cg_mismatch_is_m(self):
        # T ref, C query = generic mismatch
        arr = annotate_positions(self._pairs('T', 'C'), 'C', max_read_len=5)
        assert arr[0] == 'M'

    def test_c_to_a_mismatch_stays_m(self):
        # C ref, A query is not a bisulfite event
        arr = annotate_positions(self._pairs('C', 'A'), 'A', max_read_len=5)
        assert arr[0] == 'M'

    def test_g_to_t_mismatch_stays_m(self):
        arr = annotate_positions(self._pairs('G', 'T'), 'T', max_read_len=5)
        assert arr[0] == 'M'

    def test_no_pairs_all_dots(self):
        arr = annotate_positions([], 'ACGT', max_read_len=5)
        assert np.all(arr == '.')

    def test_pysam_lowercase_ref_handled(self):
        # pysam returns lowercase ref base for mismatches in MD tag
        pairs = [(0, 0, 'c')]
        arr = annotate_positions(pairs, 'C', max_read_len=5)
        assert arr[0] == 'C'

    def test_multiposition_sequence(self):
        # ref=CTGAC, query=TTGAC
        # pos0: C→T (conversion)
        # pos1: T→T (dot)
        # pos2: G→G (unconverted G)
        # pos3: A→A (dot)
        # pos4: C→C (unconverted C)
        ref = 'CTGAC'
        qry = 'TTGAC'
        pairs = [(i, i, r) for i, r in enumerate(ref)]
        arr = annotate_positions(pairs, qry, max_read_len=10)
        assert arr[0] == 'T'
        assert arr[1] == '.'
        assert arr[2] == 'G'
        assert arr[3] == '.'
        assert arr[4] == 'C'

    def test_output_length(self):
        arr = annotate_positions([], '', max_read_len=132)
        assert arr.shape == (132,)


# ---------------------------------------------------------------------------
# is_chimeric
# ---------------------------------------------------------------------------

class TestIsChimeric:
    def test_both_conversions_above_threshold(self):
        arr = np.full(132, '.', dtype='U1')
        arr[:3] = 'T'
        arr[3:6] = 'A'
        assert is_chimeric(arr, threshold=2) is True

    def test_only_t_conversions(self):
        arr = np.full(132, '.', dtype='U1')
        arr[:5] = 'T'
        assert is_chimeric(arr, threshold=2) is False

    def test_only_a_conversions(self):
        arr = np.full(132, '.', dtype='U1')
        arr[:5] = 'A'
        assert is_chimeric(arr, threshold=2) is False

    def test_at_threshold_not_chimeric(self):
        # exactly threshold = 2, not >, so not chimeric
        arr = np.full(132, '.', dtype='U1')
        arr[:2] = 'T'
        arr[2:4] = 'A'
        assert is_chimeric(arr, threshold=2) is False

    def test_empty_array(self):
        arr = np.full(132, '.', dtype='U1')
        assert is_chimeric(arr) is False

    def test_custom_threshold(self):
        arr = np.full(132, '.', dtype='U1')
        arr[:2] = 'T'
        arr[2:4] = 'A'
        assert is_chimeric(arr, threshold=1) is True


# ---------------------------------------------------------------------------
# minimal_unique_prefix_length
# ---------------------------------------------------------------------------

class TestMinimalUniquePrefixLength:
    def test_unique_at_length_1(self):
        assert minimal_unique_prefix_length(['ACGT', 'TGCA'], max_cycles=10) == 1

    def test_needs_length_2(self):
        assert minimal_unique_prefix_length(['ACGT', 'AGCA'], max_cycles=10) == 2

    def test_full_length_needed(self):
        assert minimal_unique_prefix_length(['AAAA', 'AAAB'], max_cycles=10) == 4

    def test_capped_by_max_cycles(self):
        # sequences differ at pos 4, but max_cycles=3 → can't distinguish
        assert minimal_unique_prefix_length(['AAAA', 'AAAB'], max_cycles=3) == 3

    def test_single_sequence(self):
        assert minimal_unique_prefix_length(['ACGT'], max_cycles=10) == 1

    def test_typical_illumina_indexes(self):
        indexes = ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA']
        length = minimal_unique_prefix_length(indexes, max_cycles=8)
        prefixes = [s[:length] for s in indexes]
        assert len(set(prefixes)) == len(indexes)

    def test_respects_sequence_length(self):
        # max_cycles longer than sequences — should not go past sequence length
        result = minimal_unique_prefix_length(['AC', 'AG'], max_cycles=100)
        assert result == 2
