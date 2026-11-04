# golooper

golooper predicts where R-loops are thermodynamically favored in a DNA sequence. It computes the equilibrium ensemble of all candidate R-loop structures — combining sequence-dependent RNA:DNA base-pairing energies with supercoiling energetics — and reports per-base probabilities, the minimum free energy structure, and ensemble properties. It is a multithreaded Go port of [rlooper2](https://github.com/garden-of-delete/rlooper2).

[![garden-of-delete](https://circleci.com/gh/garden-of-delete/go-looper.svg?style=svg)](https://circleci.com/gh/garden-of-delete/go-looper)
[![codecov](https://codecov.io/gh/garden-of-delete/go-looper/branch/main/graph/badge.svg)](https://codecov.io/gh/garden-of-delete/go-looper)

## Install

Requires Go 1.22.1 or later to compile from source. 

```
go build -o golooper        # or: make build
```

## Quick start

```
# Full ensemble analysis: per-base probability, average energy, and MFE tracks
golooper ensemble -f sequence.fa -o results --localaverageenergy

# Also call peaks and dump every structure
golooper ensemble -f sequence.fa -o results --bedfile --sensitivity=-1 --dump

# Probability of R-loop formation across a superhelicity range
golooper perloop -f sequence.fa -o results --sigma-min=-0.1 --sigma-max=-0.02

# Plasmid (circular) analysis
golooper ensemble -f plasmid.fa -o results --circular
```

## Commands

**`ensemble`** runs the full analysis on every sequence in the FASTA file and writes, using the `-o` value as a filename prefix:

| File | Contents |
|---|---|
| `_bpprob.wig` | per-base R-loop involvement probability |
| `_avgG.wig` | per-base average free energy (with `--localaverageenergy`) |
| `_mfe.wig` | 1.0 across the minimum free energy structure |
| `_bpprob.bed`, `_mfe.bed` | peak calls (with `--bedfile`) |
| `<prefix><NAME>_dump.txt` | every structure with energy and probability (with `--dump`) |

`--top N` prints the N most probable structures to stdout (residual columns are zeros unless `--residuals` is also given); `--residuals` prints ensemble-averaged residual twist, linking difference, and superhelicity.

**`perloop`** sweeps superhelical density σ over the first sequence and writes `<prefix>_perloop.txt` with one `sigma P(R-loop)` row per step. `--sigma-min` and `--sigma-max` are required; `--sigma-step` defaults to 0.01. Output flags (`--dump`, `--top`, `--bedfile`, …) are accepted but have no effect here.

## Flags

| Flag | Short | Meaning | Default |
|---|---|---|---|
| `--input` | `-f` | input FASTA file (required) | |
| `--output` | `-o` | output path prefix (required) | |
| `--sigma` | `-s` | superhelical density | −0.07 |
| `--N` | `-N` | superhelical domain size (nt), or `auto` = sequence length | 1500 (sequence length under `--circular`) |
| `--a` | `-a` | nucleation free energy, kcal/mol | 10 |
| `--minlength` | `-m` | minimum R-loop length, nt | 2 |
| `--circular` | `-C` | treat sequences as covalently closed circles | off |
| `--threads` | | worker goroutines | all CPUs |
| `--reverse` | `-r` | reverse the sequence before simulating | off |
| `--complement` | `-c` | complement the sequence before simulating | off |
| `--invert` | `-i` | reverse and complement | off |
| `--homopolymer` | `-H` | override all base-pairing energies with this value, kcal/mol | off |
| `--unconstrained` | `-u` | drop the supercoiling term (use with `--sigma 0`) | off |
| `--bedfile` | | call peaks into the bed outputs | off |
| `--sensitivity` | | peak threshold exponent: cutoff = 10^value | 1 |
| `--top` | | print the N most probable structures | 0 |
| `--dump` | `-d` | dump every structure to a file | off |
| `--residuals` | `-R` | report ensemble residual twist/linking difference/superhelicity | off |
| `--localaverageenergy` | `-l` | populate the `_avgG.wig` signal | off |

Two defaults worth knowing:

- **`--bedfile` calls no peaks at the default sensitivity.** The default threshold (10^1) exceeds any probability — preserved rlooper2 behavior. Pass e.g. `--sensitivity=-1` for a 0.1 cutoff.
- **Orientation is automatic from the header's `strand=` field** (`+` complements, `-` reverses, yielding the template strand in transcription order). `--reverse`/`--complement`/`--invert` apply on top of that.

## Input

FASTA with UCSC-style headers; only `range=` and `strand=` are parsed, and coordinates in outputs come from `range=`:

```
>NAME range=CHROM:START-END 5'pad=0 3'pad=0 strand=- repeatMasking=none
```

Only `A`/`C`/`G`/`T` (case-insensitive) are accepted.

## Circular mode

`--circular` treats each sequence as a closed circle: every rotation at every length from `--minlength` to n−1 is enumerated (full-circle structures, which would have no strand junctions, are excluded), and windows wrap through the origin. `--N` defaults to the sequence length — the plasmid is its own torsional domain — and an explicit `--N` still overrides (with a stderr note). Works identically under `perloop`.

Reading circular output: a `--dump`/`--top` row with `stop < start` is a structure wrapping through the origin. A peak spanning the origin appears as two BED rows sharing one name; pipe through `sort -k1,1 -k2,2n` before tools that require sorted BED.

## Strand separation (`golooper-sidd`)

`golooper-sidd` is a second binary, built alongside `golooper` by `make build` (or `go build -o golooper-sidd ./cmd/golooper-sidd`). It computes the stress-induced duplex destabilization (SIDD) model of Benham for a single DNA sequence: the equilibrium ensemble of locally melted regions in a superhelically stressed molecule. It reports, for each base, the probability P(x) that the base is melted and the ensemble-average free energy G(x) of opening it.

This implementation is restricted to a single denatured region per molecule: the ensemble is the ground state plus every one-bubble window (the "1R" approximation). Multiple simultaneous bubbles are not modeled.

```
# SIST-format report on stdout
golooper-sidd -f plasmid.fa --circular

# Write the report to a file and also write wig tracks
golooper-sidd -f sequence.fa --clamp -o report.txt --tracks results
```

`--tracks results` writes `results_sidd_bpprob.wig` (P(x)) and `results_sidd_gx.wig` (G(x)). The report layout matches that of a 1R-restricted SIST build (`qsidd -b -p -r` with the multi-run searches disabled): a parameter and ensemble header followed by a `Position  Base  P(x)  G(x)` table. Full SIST additionally prints multi-run counts, the total number of states, the partition function, the scaling factor, the average number of runs, the transition probability and a run time line.

SIST's recommended input range is 1,500-10,000 bp. Much longer sequences at the default sigma make every single-domain state's Boltzmann weight underflow, and the tool reports an error rather than printing wrong numbers. `--sigma` must be less than or equal to 0; SIST's positive-magnitude `-s` value corresponds to a negative `--sigma` here.

### Relationship to SIST

The calculations follow [SIST](https://github.com/CCPBioSim/SIST) (MIT license), read at commit `cc206d0`. Output is validated against a build of SIST with its two-, three- and four-run searches disabled, so that it computes the same one-run ensemble. On eight golden cases (`testdata/golden/sidd/`, regenerated by `scripts/gen-sidd-goldens.sh`, which also records the patch and build details) the output matches byte for byte.

Because the model allows one bubble per state, it is not interchangeable with full SIST when several sites compete. On pBR322 at clamp defaults, full SIST reports 79.5644 melted bases against 43.3706 from `golooper-sidd`. The site near position 307 reaches P 0.993 in full SIST and 0.891 here; the competing site near 3250 reaches 0.915 in full SIST and 0.109 here. With one bubble per state, competing sites divide the ensemble weight, so the relative P(x) between sites is not comparable to full SIST (here the weaker site is suppressed about 8x). Locations within a site agree.

That validation covers the one-run model only. It does not show that single-run results equal SIST's results. SIST sums over up to four simultaneous melted regions, and at pBR322 defaults its ensemble averages about 1.94 runs, so restricting to one run is a relevant physical approximation. Competition among multiple runs and joint competition between R-loops and bubbles are not yet implemented.

If you use `golooper-sidd`, cite the SIST paper, which its authors require for any use of these algorithms, and the paper describing the strand-separation method:

> Zhabinskaya, D., Madden, S., & Benham, C. J. (2015). SIST: stress-induced structural transitions in superhelical DNA. *Bioinformatics*, 31(3), 421-422.
>
> Fye, R. M. and Benham, C. J. (1999). Exact method for numerically analyzing a model of local denaturation in superhelically stressed DNA. *Physical Review E*, 59, 3408-3426.

(Citation strings are from the SIST README, `CCPBioSim/SIST` at `cc206d0`.)

### Flags and SIST equivalents

| `golooper-sidd` | SIST `qsidd` | Notes |
|---|---|---|
| `--infile`, `-f` | `-f` | one FASTA record; a file with no `>` header line is read as a raw sequence, as SIST does |
| `--outfile`, `-o` | stdout redirect | default is stdout |
| `--sigma` | `-s` | default −0.06. Pass a negative value for underwinding. SIST's `-s` takes a magnitude and negates positive input |
| `--temperature` | `-T` | kelvin, default 310 |
| `--salt` | `-i` | molar, default 0.01 |
| `--threshold` | `-t` | kcal/mol, default 12 |
| `--nearest-neighbor` | `-n` | default is copolymeric energetics |
| `--circular` | `-c` | |
| `--clamp` | (default geometry) | appends a 50-G clamp to a linear sequence |
| `--N` | none | superhelical domain size; see below |
| `--minlength`, `--maxlength` | none | bounds on melted-run length, defaults 1 and 250; SIST has no flag for these (its maximum window size is fixed at 250) |
| `--tracks` | none | wig output prefix |

Geometry is explicit. With neither `--clamp` nor `--circular`, the molecule is true linear DNA. SIST's default linear geometry is the clamped one, so `--clamp` is the setting that reproduces SIST's default. `--clamp` and `--circular` are mutually exclusive.

`--N auto` (the default) sets the domain size to the sequence length, or to the padded length (sequence plus clamp) under `--clamp`. An explicit `--N` is used as given, with no clamp adjustment.

Restrictions:

- Nearest-neighbor energetics require `--clamp` or `--circular`; with true-linear geometry `--nearest-neighbor` is an error.
- The input must contain exactly one sequence.
- Only `A`/`C`/`G`/`T` are accepted, as for `golooper`.

## Validation and performance

Every linear-mode output is validated against a patched build of the original C++ rlooper2 (commit `455983b`) within `max(1e-9, 1e-5·|reference|)`; goldens live in `testdata/golden/` and regenerate via `scripts/gen-golden.sh`. Circular mode is validated by independent enumeration, rotation invariance, and analytic properties. Several rlooper2 output behaviors are preserved deliberately for compatibility (the peak-threshold default and last-base signal convention among them), and a small number of behaviors depart from it where the original was internally inconsistent — exact energy ties resolve deterministically to the lowest-start window, and `--top` and `--dump` coordinates always agree.

The ensemble is O(n²) in sequence length and parallelizes across `--threads`: a 3,906 bp sequence takes ~20 ms at 8 workers (~4× over serial). `--circular` roughly doubles window count, runtime, and `--dump` size.
