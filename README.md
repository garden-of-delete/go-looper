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

## Validation and performance

Every linear-mode output is validated against a patched build of the original C++ rlooper2 (commit `455983b`) within `max(1e-9, 1e-5·|reference|)`; goldens live in `testdata/golden/` and regenerate via `scripts/gen-golden.sh`. Circular mode is validated by independent enumeration, rotation invariance, and analytic properties. Several rlooper2 output behaviors are preserved deliberately for compatibility (the peak-threshold default and last-base signal convention among them), and a small number of behaviors depart from it where the original was internally inconsistent — exact energy ties resolve deterministically to the lowest-start window, and `--top` and `--dump` coordinates always agree.

The ensemble is O(n²) in sequence length and parallelizes across `--threads`: a 3,906 bp sequence takes ~20 ms at 8 workers (~4× over serial). `--circular` roughly doubles window count, runtime, and `--dump` size.
