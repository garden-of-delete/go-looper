# golooper

golooper computes the equilibrium statistical-mechanical ensemble of R-loop
structures that a nucleic acid sequence can form. For every candidate R-loop
window (a start/end pair within the sequence), it combines a sequence-dependent
RNA:DNA base-pairing free energy with a length-dependent supercoiling free
energy (a function of superhelical density and the size of the torsionally
constrained domain) into a single free energy, then uses a Boltzmann partition
function over all windows to obtain each window's equilibrium probability, the
minimum free energy (MFE) structure, and position-resolved signals (base-pair
probability, average energy, and MFE occupancy) along the sequence. It is a Go
port of [rlooper2](https://github.com/garden-of-delete/rlooper2), aimed at predicting
where in a sequence R-loops are thermodynamically favored.

[![garden-of-delete](https://circleci.com/gh/garden-of-delete/go-looper.svg?style=svg)](https://circleci.com/gh/garden-of-delete/go-looper)
[![codecov](https://codecov.io/gh/garden-of-delete/go-looper/branch/main/graph/badge.svg)](https://codecov.io/gh/garden-of-delete/go-looper)

## Install / build

Requires Go 1.22.1 or later (the version declared in `go.mod`). From the
repository root:

```
go build -o golooper
```

or

```
make build
```

`make build` removes any previous `golooper` binary and `out/`/`test_output/`
directories, then runs `go build -o golooper`.

## Commands

### ensemble

Runs the full ensemble simulation on every sequence in a FASTA file: computes
the partition function and per-window probabilities, then writes wig/bed
tracks and optional dump, top-structure, and residual reports.

```
golooper ensemble -f <input.fa> -o <output_prefix> [flags]
```

Using `<output_prefix>` as a filename prefix, six files are always written:
`_bpprob.wig`, `_avgG.wig`, and `_mfe.wig` (UCSC `fixedStep` wiggle tracks for
base-pair probability, average energy, and MFE occupancy), `_bpprob.bed` and
`_mfe.bed` (peak calls, populated only with `--bedfile`), and `_extbpprob.wig`
(header only — this port does not implement the external/imported-structure
feature the column corresponds to). With `--dump`, an additional
`<output_prefix><sequence name>_dump.txt` is written per sequence, listing
every enumerated structure sorted by ascending energy with the ground state
first. With `--top N`, the N most probable structures are printed to stdout;
its residual columns (twist, linking difference, superhelicity) print as
zeros unless `--residuals` is also passed. With `--residuals`,
ensemble-averaged residual twist, linking difference, and superhelicity are
also printed to stdout.

### perloop

Sweeps superhelical density σ over the first sequence in a FASTA file and
reports the probability of R-loop formation at each step.

```
golooper perloop -f <input.fa> -o <output_prefix> --sigma-min <lo> --sigma-max <hi> [--sigma-step <step>] [flags]
```

`--sigma-min` and `--sigma-max` are required; `--sigma-step` defaults to
`0.01`. Writes `<output_prefix>_perloop.txt` with one row per step:
`sigma 1-P(ground state)`. `perloop` shares the full flag set with `ensemble`
but only consumes the ones that affect the model (`--sigma`, `--N`, `--a`,
`--minlength`, etc.); the ensemble-only output flags (`--dump`, `--top`,
`--residuals`, `--bedfile`, `--localaverageenergy`) are accepted but have no
effect on `perloop`'s output.

## Flags shared by both commands

| Flag | Short | Meaning |
|---|---|---|
| `--input` | `-f` | input FASTA file (required) |
| `--output` | `-o` | output path prefix (required) |
| `--sigma` | `-s` | superhelical density, e.g. `0.07` for +7% |
| `--N` | `-N` | size of the superhelicity domain in nucleotides, or `auto` to size it to the sequence length |
| `--a` | `-a` | nucleation free energy, Kcal/mol |
| `--minlength` | `-m` | minimum R-loop length, nucleotides |
| `--reverse` | `-r` | reverse the sequence (5'↔3') before simulating |
| `--complement` | `-c` | complement each base before simulating |
| `--invert` | `-i` | reverse and complement before simulating (`--reverse` + `--complement` together) |
| `--homopolymer` | `-H` | override the base-pairing energy table with this constant value, Kcal/mol |
| `--unconstrained` | `-u` | drop the supercoiling energy term from candidate structures |
| `--bedfile` | | call peaks on the wig signals and populate the bed outputs |
| `--sensitivity` | | power-of-ten exponent for the absolute peak-calling probability threshold |
| `--top` | | print the N most probable structures to stdout |
| `--dump` | `-d` | dump every enumerated structure to a file |
| `--residuals` | `-R` | compute and report residual twist/linking-difference/superhelicity |
| `--localaverageenergy` | `-l` | also compute and populate the `_avgG.wig` average-energy signal |
| `--threads` | | number of worker goroutines for the ensemble/signal computation (default: all CPUs) |
| `--circular` | `-C` | not implemented; passing it is a hard error ("circular mode is not implemented; sequences are treated as linear") |

`--reverse`/`--complement`/`--invert` apply on top of an automatic orientation
already derived from the FASTA header's `strand=` field: `+` complements the
sequence, `-` reverses it, matching the template-strand-in-transcription-order
convention described below.

`--sigma`, `--N`, `--a`, and `--minlength` are only overridden when the flag is
explicitly passed; `go run . ensemble --help` shows placeholder zero values for
these, but the model's actual defaults when a flag is omitted are σ = -0.07,
N = 1500 nt, a = 10 Kcal/mol, and minimum R-loop length = 2 nt.

`--sensitivity` sets the peak threshold to `10^sensitivity` as an absolute
probability cutoff. The default, `1`, gives a threshold of `10`, which no
probability can exceed — so by default `--bedfile` never calls any peaks. This
reproduces rlooper2's behavior deliberately; pass a negative value (e.g.
`--sensitivity -1` for a 0.1 threshold) to get peaks.

`--unconstrained` removes the supercoiling term from candidate-structure
energies but not from the ground state, matching rlooper2; it is only
physically meaningful combined with `--sigma 0`.

## FASTA header convention

Sequence headers follow the UCSC-style convention rlooper2 uses:

```
>NAME range=CHROM:START-END 5'pad=0 3'pad=0 strand=± repeatMasking=none
```

for example:

```
>HG19_AIRN_PFC53_REVERSE_dna range=PFC53FIXED:1-3908 5'pad=0 3'pad=0 strand=- repeatMasking=none
```

Only `range=` (chromosome and 1-based start/end, used for wig/bed/dump
coordinates) and `strand=` (`+` or `-`, used for the automatic orientation
above) are parsed; other fields are ignored. Only `A`/`C`/`G`/`T` bases are
accepted (case-insensitive) — anything else, including `N`, is rejected when
the file is read.

## Validation

Every energy, probability, signal, bed peak, dump row, and perloop point this
port produces is checked against a patched build of the original C++
rlooper2 at commit `455983bad1e9d49eb0e484f7746e60b16bd9868d`, at a tolerance
of `max(1e-9, 1e-5·|reference|)`. The goldens live under `testdata/golden/`
and can be regenerated at any time against a local rlooper2 checkout with
`scripts/gen-golden.sh`; the patches applied to the C++ build (documented in
`testdata/golden/README.md`) only remove undefined behavior and do not change
any computed value.

A handful of rlooper2 behaviors that affect output values are preserved
deliberately rather than "fixed": the peak-calling threshold default above,
and signals that never credit the last base of each candidate window (an
off-by-one between a window's energy length and the length used for
residual/signal bookkeeping in the original). The nucleation (junction) free
energy `--a` is likewise folded into the ground-state energy only
(`ground state = K·α²/2 − a`) and excluded from every other structure's
reported energy — so the per-structure energies in `--dump` and `--top`
output are raw supercoiling-plus-base-pairing energies, not values relative
to the ground state. There are three places this port departs from rlooper2
on purpose. First, where exact energy ties occur: rlooper2's enumeration
order resolves them arbitrarily, while golooper's parallel merge uses a
total order that always resolves a tie to the lowest-start window,
identically regardless of worker count. Second, `--top`'s coordinates on a
"-"-strand gene: rlooper2 mirrors `--top` through the sequence length while
mirroring `--dump` through the gene's genomic range, giving the two
different coordinates for the same structure; golooper uses the genomic
mirror for both, so `--top` and `--dump` always agree. Third, chromosome
names are reported as given in the FASTA header; rlooper2 lowercases the
chromosome name and then re-uppercases its last character (visible in the
golden `pfc53` case, where `PFC53FIXED` becomes `pfc53fixeD`), which this
port does not reproduce.

## Performance

Computing the ensemble is O(n²) in sequence length, since it enumerates every
`(start, end)` candidate window; runtime grows quadratically with sequence
length. Both the ensemble and signal computations parallelize over start
positions, controlled by `--threads` (default: `runtime.NumCPU()`). On a
3,906bp sequence (`res/pfc53_full.fa`), computing the ensemble takes roughly
87ms serially versus roughly 20ms at 8 workers on an Apple M2 — about a 4x
speedup — reproducible with:

```
go test ./rlooper/ -run XXX -bench BenchmarkEnsemblePfc53
```
