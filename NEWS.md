# NEWS

## v0.10.0

New opt-in options for low allele-fraction data such as cfDNA. Default output is unchanged.

- `--keep_all_ref`: always report both reference split reads (RS) and reference pairs (RP); by default one is zeroed depending on the type of alt evidence, which inflates low allele fractions (#15).
- `--clip_partner_match`: a clipped read supports the SV only if its clipped bases match the partner breakend's junction sequence better than the read's own reference. Removes clipped reads from the wrong locus (#16).
- `--clip_min_base_quality N`: ignore clips whose clipped bases have mean base quality below N (#16).
- `--fragment_counts`: new FORMAT fields counting each read pair once — `FS` split, `FP` span only, `FC` clip only (split > span > clip), `AOF`, `ROF`, and `FX` (fragments with both alt and ref evidence, counted as ref) (#17).
- `--genotype_on_fragments`: compute QA/QR, and so GT/GQ/QUAL/AB, from fragment counts; implies `--fragment_counts` (#17).

On cfDNA (off-target patients as negatives), `--clip_partner_match --fragment_counts` raises the on:off support ratio from 284 to 336 in duplex and from 69 to 75 in uncollapsed BAMs, and sensitivity at matched false-positive rate by roughly 10–20%.
