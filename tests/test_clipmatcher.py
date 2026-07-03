from .context import parsers  # ensure repo root is on sys.path for imports
import unittest
from unittest.mock import patch
import svtyper.clipmatcher as cm


class FakeFasta(object):
    """A minimal fake pysam.FastaFile that actually respects start/end
    slicing (unlike a naive stub that always returns a fixed string),
    since the anchored-junction matcher depends on correct absolute
    coordinate offsets."""

    def __init__(self, mapping):
        # mapping: {chrom: full_sequence_string}
        self.mapping = mapping

    def fetch(self, chrom, start, end):
        seq = self.mapping.get(chrom, '')
        start = max(0, start)
        end = max(start, end)
        sub = seq[start:end]
        # pad with N if the requested range runs past the end of our fake genome
        if len(sub) < end - start:
            sub = sub + 'N' * (end - start - len(sub))
        return sub

    def close(self):
        pass


def make_ref(length, inserts):
    """Build a fake reference of `length` Ns with `inserts` = {0based_start: seq} placed in."""
    chars = ['N'] * length
    for start, seq in inserts.items():
        for i, c in enumerate(seq):
            chars[start + i] = c
    return ''.join(chars)


class TestClipMatcher(unittest.TestCase):

    def test_edit_distance_max(self):
        s = 'ACGTACGT'
        t = 'ACGTTCGT'
        d = cm.edit_distance_max(s, t, 1)
        self.assertEqual(d, 1)
        d2 = cm.edit_distance_max(s, t, 0)
        self.assertTrue(d2 > 0)

    def test_effective_max_mismatch_scales_down_for_short_clips(self):
        # a flat max_mismatch of 2 is far too permissive for a 4bp clip
        # (50% mismatch tolerance); it should be scaled down to 0
        self.assertEqual(cm.effective_max_mismatch(4, 2), 0)
        # a 13bp clip should tolerate at most 1 mismatch, not the full 2
        self.assertEqual(cm.effective_max_mismatch(13, 2), 1)
        # long clips are capped at the requested max_mismatch
        self.assertEqual(cm.effective_max_mismatch(100, 2), 2)

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_match_at_anchor_forward(self, mock_fasta):
        # exact match immediately downstream of the breakpoint anchor (side A, pos=300)
        ref = make_ref(600, {300: 'AAACCCGGG'})
        mock_fasta.return_value = FakeFasta({'chr1': ref})
        clip_seq = 'AAACCCGGG'
        breakpoint = {
            'A': {'chrom': 'chr1', 'pos': 300},
            'B': {'chrom': 'chr1', 'pos': 400},
        }
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=2, max_mismatch=0)
        self.assertTrue(matched)
        self.assertEqual(dist, 0)
        self.assertEqual(side, 'A')
        self.assertEqual(orientation, 'fwd')

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_match_at_anchor_left_and_rc(self, mock_fasta):
        # exact match immediately upstream of the breakpoint anchor (side B, pos=400),
        # stored as its reverse complement in the reference so only 'rc' orientation matches
        clip_seq = 'GGTACCTA'
        upstream_seq = cm.revcomp(clip_seq)
        ref = make_ref(600, {400 - len(upstream_seq): upstream_seq})
        mock_fasta.return_value = FakeFasta({'chr1': ref})
        breakpoint = {
            'A': {'chrom': 'chr1', 'pos': 100},
            'B': {'chrom': 'chr1', 'pos': 400},
        }
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=2, max_mismatch=0)
        self.assertTrue(matched)
        self.assertEqual(dist, 0)
        self.assertEqual(side, 'B')
        self.assertEqual(orientation, 'rc')

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_incidental_match_far_from_junction_is_rejected(self, mock_fasta):
        # regression test for the "wide window" bug: an exact match of a short
        # clip placed 50bp away from the true breakpoint must NOT count, even
        # though a broad +/-100bp scan would have found it. Anchoring to a
        # small fixed radius (context) around the breakpoint's own position
        # should reject it. min_clip_length is disabled here (set to 0) so
        # this test isolates the anchoring behaviour specifically.
        clip_seq = 'CTT'
        ref = make_ref(600, {300 + 50: clip_seq})  # far from breakpoint at pos 300
        mock_fasta.return_value = FakeFasta({'chr1': ref})
        breakpoint = {
            'A': {'chrom': 'chr1', 'pos': 300},
            'B': {'chrom': 'chr1', 'pos': 400},
        }
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=5, max_mismatch=2, min_clip_length=0)
        self.assertFalse(matched)

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_short_clip_with_flat_tolerance_is_not_spuriously_matched(self, mock_fasta):
        # a 4bp clip that only shares 2/4 bases with the sequence right at the
        # anchor should NOT match once max_mismatch is scaled to clip length,
        # even though it would satisfy a flat max_mismatch=2
        ref = make_ref(600, {300: 'ATGC'})  # 2 mismatches vs 'ATCG' below
        mock_fasta.return_value = FakeFasta({'chr1': ref})
        clip_seq = 'ATCG'
        breakpoint = {
            'A': {'chrom': 'chr1', 'pos': 300},
            'B': {'chrom': 'chr1', 'pos': 400},
        }
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=2, max_mismatch=2)
        self.assertFalse(matched)

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_medium_clip_single_mismatch_within_anchor_is_found(self, mock_fasta):
        # a longer clip with a single real mismatch, located a few bp inside
        # the anchor window, should still be recovered
        ref = make_ref(600, {301: 'AAACCCGTGG'})  # true seq has G where clip has C
        mock_fasta.return_value = FakeFasta({'chr1': ref})
        clip_seq = 'AAACCCCTGG'  # single mismatch at position 7 (0-based)
        breakpoint = {
            'A': {'chrom': 'chr1', 'pos': 300},
            'B': {'chrom': 'chr1', 'pos': 400},
        }
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=2, max_mismatch=2)
        self.assertTrue(matched)
        self.assertEqual(dist, 1)

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_no_match_when_no_reference_available(self, mock_fasta):
        ref = make_ref(600, {})
        mock_fasta.return_value = FakeFasta({'chr1': ref})
        clip_seq = 'AAACCCGGG'
        breakpoint = {
            'A': {'chrom': 'chr1', 'pos': 300},
            'B': {'chrom': 'chr1', 'pos': 400},
        }
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=2, max_mismatch=0)
        self.assertFalse(matched)

    @patch('svtyper.clipmatcher.pysam.FastaFile')
    def test_short_clip_below_min_length_is_rejected_even_on_exact_match(self, mock_fasta):
        # a 2bp clip that exactly matches the reference right at the anchor
        # must still be rejected under the default min_clip_length (3), since
        # a 1-2bp exact match is essentially guaranteed to occur by chance
        # across the handful of candidate positions/directions/orientations
        # checked and carries no real information about breakpoint support.
        ref = make_ref(600, {300: 'AC'})
        mock_fasta.return_value = FakeFasta({'chr1': ref})
        clip_seq = 'AC'
        breakpoint = {
            'A': {'chrom': 'chr1', 'pos': 300},
            'B': {'chrom': 'chr1', 'pos': 400},
        }
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=2, max_mismatch=0)
        self.assertFalse(matched)

        # the same clip is allowed through (and matches) if min_clip_length is
        # explicitly lowered, confirming the rejection above was due to length
        # and not some other unrelated issue
        matched, dist, side, orientation = cm.match_clip_to_breakpoint_windows(
            clip_seq, 'fake.fa', breakpoint, context=2, max_mismatch=0, min_clip_length=2)
        self.assertTrue(matched)

    def test_ci_field_in_breakpoint_dict_is_ignored(self):
        # the matcher must not depend on breakpoint[side]['ci'] (CIPOS/CIEND)
        # being present or accurate: some SV callers don't report it, and
        # some pipelines populate it with an arbitrary placeholder value, so
        # it should have no effect on the search window. match_clip_to_breakpoint_windows
        # only requires 'chrom' and 'pos' per side.
        breakpoint_no_ci = {
            'A': {'chrom': 'chr1', 'pos': 300},
            'B': {'chrom': 'chr1', 'pos': 400},
        }
        breakpoint_with_bogus_ci = {
            'A': {'chrom': 'chr1', 'pos': 300, 'ci': [-1000, 1000]},
            'B': {'chrom': 'chr1', 'pos': 400, 'ci': [-1000, 1000]},
        }
        with patch('svtyper.clipmatcher.pysam.FastaFile') as mock_fasta:
            ref = make_ref(600, {300: 'AAACCCGGG'})
            mock_fasta.return_value = FakeFasta({'chr1': ref})
            clip_seq = 'AAACCCGGG'
            r1 = cm.match_clip_to_breakpoint_windows(clip_seq, 'fake.fa', breakpoint_no_ci, context=2, max_mismatch=0)
            r2 = cm.match_clip_to_breakpoint_windows(clip_seq, 'fake.fa', breakpoint_with_bogus_ci, context=2, max_mismatch=0)
            self.assertEqual(r1, r2)


if __name__ == '__main__':
    unittest.main()
