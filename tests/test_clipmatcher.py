from .context import parsers  # ensure repo root is on sys.path for imports
import unittest
from unittest.mock import patch
import svtyper.clipmatcher as cm


class FakeFasta(object):
    def __init__(self, mapping):
        self.mapping = mapping

    def fetch(self, chrom, start, end):
        return self.mapping.get(chrom, 'N' * 200)

    def close(self):
        pass


class TestClipMatcher(unittest.TestCase):
    def setUp(self):
        # window contains the exact sequence 'AAAACCCCTTTT' surrounded by Ns
        self.win_seq = 'N' * 50 + 'AAAACCCCTTTTGGGG' + 'N' * 50
        self.fake_fa = FakeFasta({'chr1': self.win_seq})

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_exact_match_forward(self, mock_fasta):
        mock_fasta.return_value = self.fake_fa
        clip_seq = 'AAAACCCCTTTT'
        breakpoint = {'A': {'chrom': 'chr1', 'pos': 60}, 'B': {'chrom': 'chr1', 'pos': 160}}
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=100, k=6, max_mismatch=0)
        self.assertTrue(matched)
        self.assertEqual(dist, 0)
        self.assertIn(side, ('A', 'B'))
        self.assertEqual(orientation, 'fwd')

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_rc_match(self, mock_fasta):
        mock_fasta.return_value = self.fake_fa
        clip_seq_rc = cm.revcomp('AAAACCCCTTTT')
        breakpoint = {'A': {'chrom': 'chr1', 'pos': 60}, 'B': {'chrom': 'chr1', 'pos': 160}}
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq_rc, 'fake.fa', breakpoint, context=100, k=6, max_mismatch=0)
        self.assertTrue(matched)
        self.assertEqual(dist, 0)
        self.assertIn(side, ('A', 'B'))
        self.assertEqual(orientation, 'rc')

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_no_match(self, mock_fasta):
        mock_fasta.return_value = self.fake_fa
        clip_seq = 'TTTTGGGGCCCC'  # not present in window
        breakpoint = {'A': {'chrom': 'chr1', 'pos': 60}, 'B': {'chrom': 'chr1', 'pos': 160}}
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=100, k=6, max_mismatch=0)
        self.assertFalse(matched)

    def test_edit_distance_max(self):
        s = 'ACGTACGT'
        t = 'ACGTTCGT'
        d = cm.edit_distance_max(s, t, 1)
        self.assertEqual(d, 1)
        d2 = cm.edit_distance_max(s, t, 0)
        self.assertTrue(d2 > 0)


if __name__ == '__main__':
    unittest.main()
