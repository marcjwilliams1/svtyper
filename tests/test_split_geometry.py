"""
Regression tests for SplitRead.is_split_straddle soft-clip geometry.

Focus: a soft-clip-only pseudo-split read must be recognized as supporting a
breakpoint for BOTH deletions and duplications, at BOTH breakpoint ends. DUP
junctions are everted relative to DEL, so the check against the two breakpoint
ends must be swapped for DUP soft-clips -- a regression here silently drops all
DUP clipped-read support (see is_split_straddle in svtyper/parsers.py).
"""
import unittest
import svtyper.parsers as P


class FakeRead(object):
    """Minimal stand-in for a pysam AlignedSegment, covering only the
    attributes SplitRead.is_valid() and is_split_straddle() actually read."""
    def __init__(self, reference_name, reference_start, reference_end, cigar,
                 is_reverse=False, mapq=60, qlen=80, qalen=50):
        self.query_name = "r1"
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.cigar = cigar
        self.is_reverse = is_reverse
        self.mapping_quality = mapq
        self.query_length = qlen
        self.query_alignment_length = qalen

    def has_tag(self, tag):
        return False

    def get_tag(self, tag):
        raise KeyError(tag)


SLOP = 3


def clip_supports(svtype, o1_rev, o2_rev, posA, posB, read):
    sr = P.SplitRead(read, lib=None)
    assert sr.is_valid(), "read should register as a soft-clip pseudo-split"
    assert sr.is_soft_clip, "read should be flagged is_soft_clip"
    left, right = sr.is_split_straddle("chr1", posA, [0, 0],
                                       "chr1", posB, [0, 0],
                                       o1_rev, o2_rev, svtype, SLOP)
    return bool(left or right)


class TestSplitGeometry(unittest.TestCase):
    # DEL: o1_is_reverse=False (reference_end @ posA), o2_is_reverse=True
    # (reference_start @ posB). posA/posB as received by is_split_straddle.
    def test_del_right_clip_supports_left_breakpoint(self):
        # aligned [50,100) + 30S ; aligned end == posA
        read = FakeRead("chr1", 50, 100, [(0, 50), (4, 30)])
        self.assertTrue(clip_supports("DEL", False, True, 100, 200, read))

    def test_del_left_clip_supports_right_breakpoint(self):
        # 30S + aligned [200,250) ; aligned start == posB
        read = FakeRead("chr1", 200, 250, [(4, 30), (0, 50)])
        self.assertTrue(clip_supports("DEL", False, True, 100, 200, read))

    # DUP: o1_is_reverse=True (reference_start @ posA), o2_is_reverse=False
    # (reference_end @ posB). Everted junction => swap required.
    def test_dup_right_clip_supports_right_breakpoint(self):
        # aligned [150,200) + 30S ; aligned end == posB
        read = FakeRead("chr1", 150, 200, [(0, 50), (4, 30)])
        self.assertTrue(clip_supports("DUP", True, False, 101, 200, read))

    def test_dup_left_clip_supports_left_breakpoint(self):
        # 30S + aligned [101,151) ; aligned start == posA
        read = FakeRead("chr1", 101, 151, [(4, 30), (0, 50)])
        self.assertTrue(clip_supports("DUP", True, False, 101, 200, read))

    # A clip whose aligned edge is far from either breakpoint must NOT support.
    def test_dup_clip_far_from_breakpoints_not_supported(self):
        read = FakeRead("chr1", 500, 560, [(0, 60), (4, 20)])
        self.assertFalse(clip_supports("DUP", True, False, 101, 200, read))

    # INV: an inversion breakend anchors reads by BOTH reference_end (one
    # junction) and reference_start (the other junction). o1=o2=False, posA=100,
    # posB=200. All four must be detected.
    def test_inv_reference_end_at_posA(self):
        read = FakeRead("chr1", 50, 100, [(0, 50), (4, 30)])   # reference_end=100
        self.assertTrue(clip_supports("INV", False, False, 100, 200, read))

    def test_inv_reference_end_at_posB(self):
        read = FakeRead("chr1", 150, 200, [(0, 50), (4, 30)])  # reference_end=200
        self.assertTrue(clip_supports("INV", False, False, 100, 200, read))

    def test_inv_reference_start_at_posA(self):
        read = FakeRead("chr1", 100, 150, [(4, 30), (0, 50)])  # reference_start=100
        self.assertTrue(clip_supports("INV", False, False, 100, 200, read))

    def test_inv_reference_start_at_posB(self):
        read = FakeRead("chr1", 200, 250, [(4, 30), (0, 50)])  # reference_start=200
        self.assertTrue(clip_supports("INV", False, False, 100, 200, read))

    def test_inv_clip_far_from_breakpoints_not_supported(self):
        read = FakeRead("chr1", 500, 560, [(0, 60), (4, 20)])
        self.assertFalse(clip_supports("INV", False, False, 100, 200, read))


if __name__ == "__main__":
    unittest.main()
