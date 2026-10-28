"""
Tests for the metaplot building blocks: read_utils (strand, projection, mate
merging, conversion calls) and regions (layout sniffing, parsing).

Run with:  pytest tests/
"""
import gzip

import numpy as np
import pysam
import pytest

from read_utils import conversion_calls, conversion_strand, merge_mates, project_to_window
from regions import iter_region_tuples, label_from_basename, region_layout

HEADER = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 10_000}]})


def make_read(flag=0, seq='ACGTA', cigar='5M', md='5', start=100, quals=None, xg=None):
    a = pysam.AlignedSegment(HEADER)
    a.query_name = 'r1'
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array(quals or 'I' * len(seq))
    a.set_tag('MD', md)
    if xg is not None:
        a.set_tag('XG', xg)
    return a


def arr(s):
    return np.frombuffer(s.encode('ascii'), dtype='S1').copy()


# ---------------------------------------------------------------------------
# conversion_strand
# ---------------------------------------------------------------------------

class TestConversionStrand:
    @pytest.mark.parametrize('flag, expected', [
        (99, 'OT'),   # R1 forward
        (147, 'OT'),  # R2 reverse, mate of 99
        (83, 'OB'),   # R1 reverse
        (163, 'OB'),  # R2 forward, mate of 83
        (0, 'OT'),    # single-end forward
        (16, 'OB'),   # single-end reverse
    ])
    def test_flag_rule(self, flag, expected):
        assert conversion_strand(make_read(flag=flag)) == expected

    def test_mates_agree(self):
        assert conversion_strand(make_read(flag=99)) == conversion_strand(make_read(flag=147))
        assert conversion_strand(make_read(flag=83)) == conversion_strand(make_read(flag=163))

    def test_xg_overrides_flag(self):
        # Bismark CTOT/CTOB reads: XG, not orientation, says which conversion they carry
        assert conversion_strand(make_read(flag=99, xg='GA')) == 'OB'
        assert conversion_strand(make_read(flag=83, xg='CT')) == 'OT'


# ---------------------------------------------------------------------------
# project_to_window
# ---------------------------------------------------------------------------

class TestProjectToWindow:
    def test_simple_match(self):
        bases, quals, refs = project_to_window(make_read(), win_start=98, win_len=10)
        assert bases.tobytes() == b'NNACGTANNN'
        assert refs.tobytes() == b'NNACGTANNN'
        np.testing.assert_array_equal(quals[2:7], [40] * 5)
        assert quals[:2].sum() == 0 and quals[7:].sum() == 0

    def test_ref_from_md(self):
        # read T over reference C at position 101
        bases, _, refs = project_to_window(make_read(seq='ATGTA', md='1C3'), 100, 5)
        assert bases.tobytes() == b'ATGTA'
        assert refs.tobytes() == b'ACGTA'

    def test_deletion_is_uncovered(self):
        read = make_read(seq='ACTAA', cigar='2M1D3M', md='2^G3')
        bases, quals, refs = project_to_window(read, 100, 6)
        assert bases.tobytes() == b'ACNTAA'
        assert quals[2] == 0

    def test_softclip_and_insertion_skipped(self):
        read = make_read(seq='GGACTTG', cigar='2S2M1I2M', md='4')
        bases, _, _ = project_to_window(read, 100, 4)
        assert bases.tobytes() == b'ACTG'

    def test_truncated_to_window(self):
        bases, _, _ = project_to_window(make_read(start=100), win_start=102, win_len=2)
        assert bases.tobytes() == b'GT'


# ---------------------------------------------------------------------------
# merge_mates (amplicon-smf fork rules)
# ---------------------------------------------------------------------------

class TestMergeMates:
    def mate(self, seq, quals, ref=None):
        return arr(seq), np.asarray(quals, dtype=int), arr(ref or seq)

    def test_rules(self):
        #        pos:     0    1    2    3    4    5
        m1 = self.mate('NCCTAG', [0, 0, 30, 30, 30, 30], ref='NCCCAG')
        m2 = self.mate('NNTCAA', [0, 0, 20, 30, 30, 30], ref='NNCCAG')
        bases, refs = merge_mates(m1, m2)
        assert bases[0] == b'N'   # 1: uncovered
        assert bases[1] == b'N'   # 2: one mate, quality 0
        assert bases[2] == b'C'   # 4: different qualities -> higher (m1, q30)
        assert bases[3] == b'N'   # 6: same quality, different bases
        assert bases[4] == b'A'   # 5: same quality, same base
        assert bases[5] == b'N'   # 6: G vs A at equal quality
        assert refs.tobytes() == b'NCCCAG'

    def test_one_mate_only(self):
        m1 = self.mate('NNCT', [0, 0, 30, 30])
        m2 = self.mate('NNNN', [0, 0, 0, 0])
        bases, _ = merge_mates(m1, m2)
        assert bases.tobytes() == b'NNCT'   # 3: single coverage, q > 0

    def test_single_end(self):
        bases, refs = merge_mates(self.mate('ACGT', [30, 0, 30, 30]))
        assert bases.tobytes() == b'ANGT'
        assert refs.tobytes() == b'ACGT'

    def test_overlap_counted_once(self):
        # both mates cover the same C->T; the merged molecule has it once
        m = self.mate('TTTT', [30] * 4, ref='CCCC')
        bases, refs = merge_mates(m, m)
        covered, converted = conversion_calls('OT', refs, bases)
        assert covered.sum() == 4 and converted.sum() == 4


# ---------------------------------------------------------------------------
# conversion_calls
# ---------------------------------------------------------------------------

class TestConversionCalls:
    def test_ot(self):
        covered, converted = conversion_calls('OT', arr('CCCGA'), arr('CTAGA'))
        np.testing.assert_array_equal(covered, [True, True, False, False, False])
        np.testing.assert_array_equal(converted, [False, True, False, False, False])

    def test_ob(self):
        covered, converted = conversion_calls('OB', arr('GGGCA'), arr('GACCA'))
        np.testing.assert_array_equal(covered, [True, True, False, False, False])
        np.testing.assert_array_equal(converted, [False, True, False, False, False])

    def test_n_not_informative(self):
        covered, _ = conversion_calls('OT', arr('C'), arr('N'))
        assert not covered.any()


# ---------------------------------------------------------------------------
# regions
# ---------------------------------------------------------------------------

def write(path, lines, gz=False):
    opener = gzip.open if gz else open
    with opener(path, 'wt') as fh:
        fh.write('\n'.join(lines) + '\n')
    return str(path)


class TestRegions:
    def test_bed6_strand_is_col5(self, tmp_path):
        p = write(tmp_path / 'sites.bed', ['chr1\t10\t20\tname\t0\t-'])
        assert region_layout(p) == ((0, 1, 2, 5), 'sites')
        assert list(iter_region_tuples(p)) == [('chr1', '10', '20', '-')]

    def test_4col(self, tmp_path):
        p = write(tmp_path / 'a.bed', ['#header', 'chr1\t5\t5\t+'])
        assert list(iter_region_tuples(p)) == [('chr1', '5', '5', '+')]

    def test_refflat_tss_7col(self, tmp_path):
        p = write(tmp_path / 'TSS.bed', ['#chr\tleft\tright\tstrand\tg\tid\ttx',
                                         'chr1\t100\t100\t-\tGENE\tID\tTX'])
        assert region_layout(p)[0] == (0, 1, 2, 3)

    def test_11col_gzip_by_magic(self, tmp_path):
        fields = ['site', 'CTCF', 'chr1', '100', '119', '+'] + ['x'] * 5
        p = write(tmp_path / 'motifs.txt', ['\t'.join(fields)], gz=True)  # no .gz suffix
        assert region_layout(p) == ((2, 3, 4, 5), 'CTCF')

    def test_22col(self, tmp_path):
        fields = ['p'] * 11 + ['site', 'GATA1', 'chr2', '50', '55', '-'] + ['x'] * 5
        p = write(tmp_path / 'peaks.100radius', ['\t'.join(fields)])
        assert list(iter_region_tuples(p)) == [('chr2', '50', '55', '-')]

    def test_unknown_layout_raises(self, tmp_path):
        p = write(tmp_path / 'bad.bed', ['\t'.join(['x'] * 16)])
        with pytest.raises(ValueError, match='unrecognised column count 16'):
            region_layout(p)

    def test_label_from_basename(self):
        assert label_from_basename('/a/TSS_K562_RPKM.Q1.bed') == 'TSS_K562_RPKM.Q1'
        assert label_from_basename('/a/sites.bed.gz') == 'sites'
