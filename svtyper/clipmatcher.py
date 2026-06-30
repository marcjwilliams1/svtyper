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


def _kmer_prefilter(clip_seq, window_seq, k=8):
    """Return True if any k-mer from clip_seq appears in window_seq."""
    if clip_seq is None or window_seq is None:
        return False
    L = len(clip_seq)
    k = min(k, max(1, L))
    clip_seq = clip_seq.upper()
    window_seq = window_seq.upper()
    seen = set()
    for i in range(0, L - k + 1):
        seen.add(clip_seq[i:i+k])
    for km in seen:
        if km in window_seq:
            return True
    return False


def match_clip_to_breakpoint_windows(clip_seq, ref_fasta_path, breakpoint, context=100, k=8, max_mismatch=2):
    """
    Attempt to match clip_seq anywhere in +/-context bases around each breakpoint side.
    Returns (matched, distance, side, orientation) where side is 'A' or 'B' and orientation is 'rc' or 'fwd'.
    Conservative: returns False if reference not available or no match.
    """
    if clip_seq is None or len(clip_seq) == 0:
        return (False, None, None, None)

    clip_seq_u = clip_seq.upper()
    clip_rc = revcomp(clip_seq_u)

    # open fasta
    try:
        fasta = pysam.FastaFile(ref_fasta_path)
    except Exception:
        return (False, None, None, None)

    L = len(clip_seq_u)
    maxd = max_mismatch

    for side in ('A', 'B'):
        chrom = breakpoint[side]['chrom']
        pos = int(breakpoint[side]['pos'])
        # fetch a window around pos
        start = max(0, pos - context)
        end = pos + context
        try:
            win = fasta.fetch(chrom, start, end).upper()
        except Exception:
            continue

        # quick k-mer prefilter for forward and rc
        if not _kmer_prefilter(clip_seq_u, win, k) and not _kmer_prefilter(clip_rc, win, k):
            continue

        # sliding window exact-length alignment with edit-distance cutoff
        if len(win) < L:
            # compare whole window to clip_seq prefix/suffix
            d1 = edit_distance_max(clip_seq_u[:len(win)], win, maxd)
            if d1 <= maxd:
                fasta.close()
                return (True, d1, side, 'fwd')
            d2 = edit_distance_max(clip_rc[:len(win)], win, maxd)
            if d2 <= maxd:
                fasta.close()
                return (True, d2, side, 'rc')
            continue

        # iterate candidate substrings
        for i in range(0, len(win) - L + 1):
            sub = win[i:i+L]
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
