from .context import classic  # ensure repo root is on sys.path for imports
import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import svtyper.clipmatcher as cm
from .test_clipmatcher import FakeFasta

HERE = os.path.dirname(__file__)
FIX = os.path.join(HERE, "data", "clipmatch")

# synthetic genome: two contigs of distinct, non-repetitive sequence
A_SEQ = ("ACGTTGCAAGCTTGACCTAGGATCCGATTACAGGCATCGTAGCTAGCTTACGGATCCAATGCGTACG" * 4)[:240]
B_SEQ = ("TTGGCCAATCGGATAGCTCGAGGTACCTTAGGCAATCTCGACGTCATGGAACTTCGGTACAGTCAAG" * 4)[:240]


def bp(a_rev, b_rev, svtype="BND", pos_a=100, pos_b=150):
    return {"A": {"chrom": "a", "pos": pos_a, "ci": [0, 0], "is_reverse": a_rev},
            "B": {"chrom": "b", "pos": pos_b, "ci": [0, 0], "is_reverse": b_rev},
            "svtype": svtype}


class TestMatchClipToPartner(unittest.TestCase):
    """A forward breakend at a:100 (sequence kept to its left) joined to a reverse
    breakend at b:150 (sequence kept from 150 on): the junction reads
    a[..100] + b[150..], so a read ending at a:100 clips with b's sequence."""

    def setUp(self):
        self.fasta = FakeFasta({"a": A_SEQ, "b": B_SEQ})

    def match(self, clip, side, edge, breakpoint, **kw):
        with patch("svtyper.clipmatcher.pysam.FastaFile", return_value=self.fasta):
            return cm.match_clip_to_partner(clip, side, "a", edge, "ref.fa", breakpoint, **kw)

    def test_partner_clip_matches(self):
        clip = B_SEQ[149:164]                              # what the junction adds
        ok, d_partner, d_own, side = self.match(clip, "right", 100, bp(False, True))
        self.assertTrue(ok)
        self.assertEqual((d_partner, side), (0, "A"))

    def test_partner_clip_with_mismatch_matches(self):
        clip = B_SEQ[149:156] + ("A" if B_SEQ[156] != "A" else "C") + B_SEQ[157:164]
        self.assertTrue(self.match(clip, "right", 100, bp(False, True))[0])

    def test_own_reference_clip_rejected(self):
        # a read clipped for another reason: its clip is the reference just past the edge
        clip = A_SEQ[100:115]
        ok, d_partner, d_own, _ = self.match(clip, "right", 100, bp(False, True))
        self.assertFalse(ok)
        self.assertEqual(d_own, 0)

    def test_old_matcher_accepts_own_reference_clip(self):
        # documents the behaviour --clip_partner_match exists to fix
        with patch("svtyper.clipmatcher.pysam.FastaFile", return_value=self.fasta):
            matched = cm.match_clip_to_breakpoint_windows(A_SEQ[100:115], "ref.fa", bp(False, True))[0]
        self.assertTrue(matched)

    def test_wrong_side_rejected(self):
        # a forward breakend's junction is to its right: a left clip can't support it
        clip = B_SEQ[149:164]
        self.assertFalse(self.match(clip, "left", 100, bp(False, True))[0])

    def test_not_near_a_breakend_rejected(self):
        self.assertFalse(self.match(B_SEQ[149:164], "right", 130, bp(False, True))[0])

    def test_forward_partner_is_reverse_complemented(self):
        # partner kept up to b:150 and joined inverted: the junction adds revcomp(b[..150])
        clip = cm.revcomp(B_SEQ[135:150])
        self.assertTrue(self.match(clip, "right", 100, bp(False, False))[0])
        self.assertFalse(self.match(B_SEQ[149:164], "right", 100, bp(False, False))[0])

    def test_reverse_own_breakend_left_clip(self):
        # own breakend reverse (kept from a:100 on), partner forward: the junction reads
        # b[..150] + a[100..], so a read starting at a:100 left-clips with b[..150]
        clip = B_SEQ[135:150]
        self.assertTrue(self.match(clip, "left", 99, bp(True, False))[0])

    def test_homology_tie_rejected(self):
        # partner sequence identical to the own continuation: the clip carries no information
        b = B_SEQ[:149] + A_SEQ[100:140] + B_SEQ[189:]
        fasta = FakeFasta({"a": A_SEQ, "b": b})
        with patch("svtyper.clipmatcher.pysam.FastaFile", return_value=fasta):
            ok = cm.match_clip_to_partner(A_SEQ[100:115], "right", "a", 100, "ref.fa", bp(False, True))[0]
        self.assertFalse(ok)

    def test_inv_allows_both_sides(self):
        # INV ++ junction: a[..100] + revcomp(b[..150]); -- junction: revcomp(a[100..]) + b[150..]
        self.assertTrue(self.match(cm.revcomp(B_SEQ[135:150]), "right", 100, bp(False, False, "INV"))[0])
        self.assertTrue(self.match(cm.revcomp(B_SEQ[149:164]), "left", 99, bp(False, False, "INV"))[0])

    def test_short_clip_rejected(self):
        self.assertFalse(self.match(B_SEQ[149:156], "right", 100, bp(False, True))[0])


class TestClipQuality(unittest.TestCase):

    def read(self, cigar, quals):
        return SimpleNamespace(cigartuples=cigar, query_qualities=quals)

    def test_mean_quality_threshold(self):
        r = self.read([(0, 5), (4, 3)], [30] * 5 + [10, 12, 14])
        self.assertFalse(cm.clip_quality_ok(r, "right", 20))
        self.assertTrue(cm.clip_quality_ok(r, "right", 12))

    def test_left_clip(self):
        r = self.read([(4, 2), (0, 5)], [35, 37] + [10] * 5)
        self.assertTrue(cm.clip_quality_ok(r, "left", 30))

    def test_nothing_to_judge_passes(self):
        self.assertTrue(cm.clip_quality_ok(self.read([(5, 4), (0, 5)], [10] * 5), "left", 30))   # hard clip
        self.assertTrue(cm.clip_quality_ok(self.read([(0, 5)], None), "right", 30))


class TestClipOptionsOnFixture(unittest.TestCase):
    """svtyper on tests/data/clipmatch (see its README): sv1 a real junction, sv2 clips
    that are just the reads' own reference, sv3 low-quality clips."""

    def counts(self, **opts):
        with tempfile.TemporaryDirectory() as tmp:
            out_vcf = os.path.join(tmp, "out.vcf")
            with open(os.path.join(FIX, "svs.vcf")) as inf, open(out_vcf, "w") as outf:
                classic.sv_genotype(bam_string=os.path.join(FIX, "reads.bam"), vcf_in=inf, vcf_out=outf,
                                    min_aligned=25, split_weight=1, disc_weight=1, num_samp=1000000,
                                    lib_info_path=os.path.join(FIX, "lib.json"), debug=False,
                                    clip_read_support=True, ref_fasta=os.path.join(FIX, "ref.fa"),
                                    both_sides=True, keep_all_ref=True, **opts)
            lines = open(out_vcf).read().splitlines()
        res = {}
        for line in lines:
            if line.startswith("#"):
                continue
            t = line.split("\t")
            fmt = dict(zip(t[8].split(":"), t[-1].split(":")))
            res[t[2].rsplit("_", 1)[0]] = {k: int(fmt[k]) for k in ("AS", "ASC", "AP")}
        return res

    def test_default_matcher(self):
        c = self.counts()
        self.assertEqual(c["sv1"], {"AS": 8, "ASC": 6, "AP": 11})
        self.assertGreater(c["sv2"]["ASC"], 0)
        self.assertGreater(c["sv3"]["ASC"], 0)

    def test_partner_match_keeps_junction_drops_own_reference_clips(self):
        c = self.counts(clip_partner_match=True)
        self.assertEqual(c["sv1"], {"AS": 8, "ASC": 6, "AP": 11})
        self.assertEqual(c["sv2"]["ASC"], 0)

    def test_min_base_quality_drops_low_quality_clips(self):
        c = self.counts(clip_min_base_quality=30)
        self.assertEqual(c["sv1"], {"AS": 8, "ASC": 6, "AP": 11})
        self.assertEqual(c["sv3"]["ASC"], 0)
        self.assertGreater(c["sv2"]["ASC"], 0)


if __name__ == "__main__":
    unittest.main()
