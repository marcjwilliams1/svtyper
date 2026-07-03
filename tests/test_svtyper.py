
from .context import parsers as p
from .context import classic
import unittest, os, subprocess

HERE = os.path.dirname(__file__)
in_vcf = os.path.join(HERE, "data/example.vcf")
in_bam = os.path.join(HERE, "data/NA12878.target_loci.sorted.bam")
lib_info_json = os.path.join(HERE, "data/NA12878.bam.json")
out_vcf = os.path.join(HERE, "data/out.vcf")
expected_out_vcf = os.path.join(HERE, "data/example.gt.vcf")

class TestCigarParsing(unittest.TestCase):
    def test_cigarstring_to_tuple(self):
        string1 = '5H3S2D1N5M3I2P2X1='
        self.assertEqual(p.SplitRead.cigarstring_to_tuple(string1),
                [(5, 5), (4, 3), (2, 2), (3, 1),
                    (0, 5), (1, 3), (6, 2), (8, 2),
                    (7, 1)])

    def test_get_query_pos_from_cigar(self):
        # forward
        cigar_string = '2S3M1D2M2I3M3S'
        cigar = p.SplitRead.cigarstring_to_tuple(cigar_string)
        query_pos = p.SplitRead.SplitPiece.get_query_pos_from_cigar(cigar, True)
        self.assertEqual(query_pos.query_start, 3)
        self.assertEqual(query_pos.query_end, 13)
        self.assertEqual(query_pos.query_length, 15)

        # get_query_pos_from_cigar currently modifies the cigar list in place.
        # that's why the code below doesn't work as intended.
        query_pos = p.SplitRead.SplitPiece.get_query_pos_from_cigar(cigar, False)
        self.assertEqual(query_pos.query_start, 2)
        self.assertEqual(query_pos.query_end, 12)
        self.assertEqual(query_pos.query_length, 15)

    def test_get_reference_end_from_cigar(self):
        cigar_string = '2S5M3D2M3S'
        cigar = p.SplitRead.cigarstring_to_tuple(cigar_string)
        self.assertEqual(p.SplitRead.get_reference_end_from_cigar(1, cigar), 11)

    def test_get_start_diagonal(self):
        cigar_string = '2S5M3D1I1M3S'
        split_piece = p.SplitRead.SplitPiece(1, 25, True, p.SplitRead.cigarstring_to_tuple(cigar_string), 60)
        self.assertEqual(p.SplitRead.get_start_diagonal(split_piece), 23)
        split_piece2 = p.SplitRead.SplitPiece(1, 25, False, p.SplitRead.cigarstring_to_tuple(cigar_string), 60)
        self.assertEqual(p.SplitRead.get_start_diagonal(split_piece2), 23)

    def test_get_end_diagonal(self):
        cigar_string = '2S5M3D2I1M3S'
        split_piece = p.SplitRead.SplitPiece(1, 25, True, p.SplitRead.cigarstring_to_tuple(cigar_string), 60)
        split_piece.set_reference_end(34)
        self.assertEqual(p.SplitRead.get_end_diagonal(split_piece), 34 - (2 + 8))
        split_piece2 = p.SplitRead.SplitPiece(1, 25, False, p.SplitRead.cigarstring_to_tuple(cigar_string), 60)
        split_piece2.set_reference_end(34)
        self.assertEqual(p.SplitRead.get_end_diagonal(split_piece2), 34 - (2 + 8))

_REF_CONSUMING_OPS = (0, 2, 3, 7, 8)  # M, D, N, =, X


class FakeRead(object):
    """Minimal stand-in for a pysam.AlignedSegment, exposing just the
    attributes get_clipped_sequence() reads."""
    def __init__(self, cigar_string, seq, reference_start=1000, reference_name='1'):
        self.cigar = p.SplitRead.cigarstring_to_tuple(cigar_string)
        self.query_sequence = seq
        self.query_name = 'fake_read'
        self.reference_start = reference_start
        self.reference_name = reference_name
        ref_span = sum(n for op, n in self.cigar if op in _REF_CONSUMING_OPS)
        self.reference_end = reference_start + ref_span


class TestGetClippedSequence(unittest.TestCase):
    def _clip(self, cigar_string, seq, reference_start=1000, reference_name='1', anchor_positions=None):
        split = p.SplitRead(FakeRead(cigar_string, seq, reference_start, reference_name), None)
        return split.get_clipped_sequence(anchor_positions=anchor_positions)

    def test_left_clip_only(self):
        self.assertEqual(self._clip('3S7M', 'AAATTTTTTT'), ('AAA', 'left'))

    def test_right_clip_only(self):
        self.assertEqual(self._clip('7M3S', 'TTTTTTTAAA'), ('AAA', 'right'))

    def test_no_clip(self):
        self.assertEqual(self._clip('10M', 'TTTTTTTTTT'), (None, None))

    def test_both_clipped_no_anchor_prefers_longer_right(self):
        # left clip is 3bp, right clip is 5bp: without anchor info, fall
        # back to the longer, more informative clip (right) rather than
        # always the left one
        self.assertEqual(self._clip('3S7M5S', 'AAATTTTTTTGGGGG'), ('GGGGG', 'right'))

    def test_both_clipped_no_anchor_prefers_longer_left(self):
        self.assertEqual(self._clip('5S7M3S', 'GGGGGTTTTTTTAAA'), ('GGGGG', 'left'))

    def test_both_clipped_no_anchor_equal_length_ties_to_left(self):
        self.assertEqual(self._clip('3S7M3S', 'AAATTTTTTTGGG'), ('AAA', 'left'))

    def test_both_clipped_prefers_side_closest_to_anchor(self):
        # reference_start=1000, reference_end=1007 (3S7M5S).
        # left edge (1000) is right next to breakpoint anchor 999: pick left,
        # even though the right clip (5bp) is longer than the left (3bp).
        self.assertEqual(
            self._clip('3S7M5S', 'AAATTTTTTTGGGGG', reference_start=1000, anchor_positions=[('1', 999)]),
            ('AAA', 'left'))

    def test_both_clipped_prefers_side_closest_to_anchor_right(self):
        # right edge (1007) is right next to breakpoint anchor 1008: pick
        # right, even though it's not the longer clip in this case either
        # (both 3bp here, but distance should still decide it explicitly).
        self.assertEqual(
            self._clip('3S7M3S', 'AAATTTTTTTGGG', reference_start=1000, anchor_positions=[('1', 1008)]),
            ('GGG', 'right'))

    def test_anchor_on_different_chrom_is_ignored(self):
        # anchor is on chrom '2', read is on chrom '1' -> no usable anchor,
        # falls back to longer-clip behavior (right, 5bp > 3bp)
        self.assertEqual(
            self._clip('3S7M5S', 'AAATTTTTTTGGGGG', reference_start=1000, reference_name='1',
                        anchor_positions=[('2', 999)]),
            ('GGGGG', 'right'))

    def test_anchor_equidistant_falls_back_to_longer_clip(self):
        # anchor is exactly between the two edges (1000 and 1007) -> tie,
        # falls back to the longer clip (right, 5bp)
        self.assertEqual(
            self._clip('3S7M5S', 'AAATTTTTTTGGGGG', reference_start=1000, anchor_positions=[('1', 1003.5)]),
            ('GGGGG', 'right'))


class TestIntegration(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        if os.path.exists(out_vcf):
            os.remove(out_vcf)

    def test_integration(self):
        with open(in_vcf, "r") as inf, open(out_vcf, "w") as outf:
            classic.sv_genotype(bam_string=in_bam,
                             vcf_in=inf,
                             vcf_out=outf,
                             min_aligned=20,
                             split_weight=1,
                             disc_weight=1,
                             num_samp=1000000,
                             lib_info_path=lib_info_json,
                             debug=False,
                             alignment_outpath=None,
                             ref_fasta=None,
                             sum_quals=False,
                             max_reads=None,
                             max_ci_dist=1e10)

        fail_msg = "did not file output vcf '{}' after running sv_genotype".format(out_vcf)
        self.assertTrue(os.path.exists(out_vcf), fail_msg)

        fail_msg = ("output vcf '{}' "
                    "did not match expected "
                    "output vcf '{}'").format(out_vcf, expected_out_vcf)
        self.assertTrue(self.diff(), fail_msg)

    def diff(self):
        cmd = ['diff', "-I", "^##fileDate=", expected_out_vcf, out_vcf]

        rv = None
        with open(os.devnull, "w") as f:
            rv = subprocess.call(cmd, stdout=f)

        result = rv == 0
        return result

if __name__ == '__main__':
    unittest.main(verbosity=2)
