from __future__ import print_function

import pysam

# simple helpers for clipped-read matching

_base_comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def revcomp(seq):
    try:
        return seq.translate(_base_comp)[::-1]
    except Exception:
        # fallback for Python2-style environments
        trans = {}
        for a,b in zip('ACGTNacgtn','TGCANtgcan'):
            trans[ord(a)] = ord(b)
        return seq.translate(trans)[::-1]


def edit_distance_max(s, t, maxd):
    """
    Compute Levenshtein distance with early bail-out if distance exceeds maxd.
    Returns an integer > maxd if exceeded.
    """
    if s is None or t is None:
        return maxd + 1
    s = s.upper()
    t = t.upper()
    ls = len(s)
    lt = len(t)
    # quick bound
    if abs(ls - lt) > maxd:
        return maxd + 1

    # ensure s is the shorter
    if ls > lt:
        s, t = t, s
        ls, lt = lt, ls

    previous = list(range(ls + 1))
    for j in range(1, lt + 1):
        c = t[j - 1]
        current = [j] + [0] * ls
        # we can track the minimal value in this row to early abort
        row_min = current[0]
        for i in range(1, ls + 1):
            add = 1 if s[i - 1] != c else 0
            current[i] = min(previous[i] + 1,      # deletion
                             current[i - 1] + 1,   # insertion
                             previous[i - 1] + add) # substitution
            if current[i] < row_min:
                row_min = current[i]
        if row_min > maxd:
            return maxd + 1
        previous = current
    dist = previous[-1]
    return dist


def effective_max_mismatch(clip_len, max_mismatch):
    """
    Scale down the allowed edit distance for short clips so that the tolerance
    never exceeds a fixed fraction of the clip length. Without this, a flat
    max_mismatch (e.g. 2) is nearly meaningless for very short clips (e.g. a
    4bp clip with 2 allowed mismatches is only 50% identity) and produces
    spurious matches almost anywhere in the search window.
    """
    if clip_len <= 0:
        return 0
    # allow roughly 1 mismatch per 7 bases, capped by the requested max_mismatch
    scaled = clip_len // 7
    return max(0, min(max_mismatch, scaled))


def match_clip_to_breakpoint_windows(clip_seq, ref_fasta_path, breakpoint, context=5, k=8, max_mismatch=2, min_clip_length=3):
    """
    Match a clipped read's sequence against the reference immediately
    adjacent to the breakpoint junction, rather than scanning broadly across
    a wide window around the breakpoint.

    Rationale: searching a broad window (e.g. +/-100bp) for *any* occurrence
    of the clip sequence is not the same as testing whether the clip
    actually represents the breakpoint-adjacent sequence. Short clips in
    particular have few possible sequences (e.g. only 4^3=64 possible 3-mers)
    so an incidental exact match somewhere unrelated within a wide window is
    expected by chance alone, even though it says nothing about whether the
    read truly supports the junction. (Confirmed empirically: some clips
    previously "matched" 50-90bp away from the true breakpoint.)

    Instead, for each breakpoint side we only consider anchor positions
    within `pos +/- context` (default 5bp) of the breakpoint's own reported
    position. This intentionally does NOT use the VCF's CIPOS/CIEND
    confidence interval: not all SV callers report it reliably (or at all),
    and some pipelines pad it with an arbitrary/placeholder value, so it
    can't be trusted as a meaningful bound on positional uncertainty.
    `context` alone defines the full search radius around `pos`.

    For each anchor position we check the reference immediately to its left
    and right (since a clip can represent sequence that continues in either
    direction relative to the anchor), in both forward and
    reverse-complement orientation (since the correct strand/direction isn't
    always known without deeper breakend bookkeeping).

    Clips shorter than `min_clip_length` are rejected outright, regardless of
    match quality: very short clips (1-2bp) have so few possible sequences
    (4 or 16) that even an *exact* match within a small anchor window is
    expected to occur by chance alone across the handful of candidate
    positions/directions/orientations checked, so they carry essentially no
    information about whether the read truly supports the breakpoint.

    `k` is accepted for backward compatibility with existing CLI wiring but
    is no longer used: the k-mer prefilter is unnecessary now that the
    candidate search space is already small.

    Returns (matched, distance, side, orientation).
    """
    if clip_seq is None or len(clip_seq) == 0:
        return (False, None, None, None)
    if len(clip_seq) < min_clip_length:
        return (False, None, None, None)

    clip_seq_u = clip_seq.upper()
    clip_rc = revcomp(clip_seq_u)
    L = len(clip_seq_u)
    maxd = effective_max_mismatch(L, max_mismatch)

    try:
        fasta = pysam.FastaFile(ref_fasta_path)
    except Exception:
        return (False, None, None, None)

    for side in ('A', 'B'):
        info = breakpoint.get(side)
        if not info:
            continue
        chrom = info['chrom']
        pos = int(info['pos'])

        # candidate anchor positions: a fixed radius (`context`) around the
        # breakpoint's own position -- deliberately independent of CIPOS/CIEND
        anchor_lo = pos - context
        anchor_hi = pos + context

        # fetch one block covering all candidate positions (and enough
        # padding on either side for the clip length itself)
        block_start0 = max(0, anchor_lo - L - 1)  # 0-based
        block_end0 = anchor_hi + L                # 0-based, exclusive
        try:
            block = fasta.fetch(chrom, block_start0, block_end0).upper()
        except Exception:
            continue

        for apos in range(anchor_lo, anchor_hi + 1):
            # "right": reference bases immediately downstream of apos
            r_off = apos - block_start0
            sub_right = block[r_off:r_off + L]
            # "left": reference bases immediately upstream of (and including) apos
            l_off = apos - L - block_start0
            sub_left = block[max(0, l_off):max(0, l_off) + L] if l_off >= 0 else None

            for sub in (sub_right, sub_left):
                if not sub or len(sub) != L:
                    continue
                d = edit_distance_max(clip_seq_u, sub, maxd)
                if d <= maxd:
                    fasta.close()
                    return (True, d, side, 'fwd')
                d = edit_distance_max(clip_rc, sub, maxd)
                if d <= maxd:
                    fasta.close()
                    return (True, d, side, 'rc')

    fasta.close()
    return (False, None, None, None)
