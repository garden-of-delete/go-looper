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
| `--N` | `-N` | size of the superhelicity domain in nucleotides, or `auto` to size it to the sequence length (the default under `--circular`) |
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
| `--circular` | `-C` | treat the sequence as a covalently closed circle (see Circular mode below) |

`--reverse`/`--complement`/`--invert` apply on top of an automatic orientation
already derived from the FASTA header's `strand=` field: `+` complements the
sequence, `-` reverses it, matching the template-strand-in-transcription-order
convention described below.

`--sigma`, `--N`, `--a`, and `--minlength` are only overridden when the flag is
explicitly passed; `go run . ensemble --help` shows placeholder zero values for
these, but the model's actual defaults when a flag is omitted are σ = -0.07,
N = 1500 nt (or the sequence length under `--circular`), a = 10 Kcal/mol, and
minimum R-loop length = 2 nt.

`--sensitivity` sets the peak threshold to `10^sensitivity` as an absolute
probability cutoff. The default, `1`, gives a threshold of `10`, which no
probability can exceed — so by default `--bedfile` never calls any peaks. This
reproduces rlooper2's behavior deliberately; pass a negative value (e.g.
`--sensitivity -1` for a 0.1 threshold) to get peaks.

`--unconstrained` removes the supercoiling term from candidate-structure
energies but not from the ground state, matching rlooper2; it is only
physically meaningful combined with `--sigma 0`.

## Circular mode

`--circular` treats a sequence as a covalently closed circle of length n
instead of a linear molecule, for both `ensemble` and `perloop`. The
structure ensemble becomes every window (start, length) with start in
`[0, n-1]` and length in `[minlength, n-1]` — every rotation at every length
is enumerated, for `n·(n-minlength)` structures in total. The full-circle
window (length n) is excluded: every R-loop structure has two real strand
junctions, but a window spanning the entire circle would have none, so that
length is never considered. A window whose start+length exceeds n wraps
through the origin and credits bases `(start+j) mod n` for `j` in
`[0, length-2]` — the same "last base of the window is uncredited" convention
linear mode uses, applied cyclically.

This is Go-defined semantics, not a port of the C++ rlooper2's `--circular`.
That implementation hangs for `--minlength >= 3` and omits every wrapped
structure from its per-base signals while still counting those structures in
the partition function — it is not a usable reference and is not treated as
one here. See
`docs/superpowers/specs/2026-11-03-circular-mode-design.md` for the full
semantics and the empirical trace of the C++ behavior. Circular-mode output
is validated independently (below) rather than against a C++ golden.

Under `--circular`, `--N` defaults per gene to the sequence length rather
than the model's usual 1500nt default — the closed circle is its own
torsional domain — using the same path as `--N auto`. An explicit numeric
`--N` overrides this default; if it differs from the sequence length, a note
is printed to stderr: `note: --N %d differs from circular sequence length
%d`. `--N auto` is also accepted under `--circular`, redundantly setting the
same default.

`--dump` and `--top` rows for a wrapped structure print with
`stop_position < start_position` in genomic coordinates; this is how the
wrap is encoded. A row with `stop_position >= start_position` is a
non-wrapped structure, as in linear mode. `--bedfile` peak-calling treats
the signal as cyclic: if the runs touching index 0 and index n-1 are both
above threshold, they are one peak spanning the origin, emitted as two BED
rows — the arc ending at n-1, then the arc starting at 0 — sharing one peak
name/index, so the two rows can be rejoined downstream by name. If the
entire signal is above threshold, it is emitted as a single full-span row.
Because the end-arc row is emitted before the start-arc row by design, an
origin-spanning peak's two rows are not in ascending genomic order; pipe
circular `--bedfile` output through `sort -k1,1 -k2,2n` before feeding it to
tools that assume BED-sorted input.

`perloop` accepts `--circular` identically to `ensemble`. A side effect:
`--N auto` (and the circular `--N` default) now takes effect under
`perloop`, in both circular and linear runs. In the C++ rlooper2,
`simulation_B` (perloop's equivalent) silently ignored `auto_domain_size`;
that was untested behavior rather than a deliberate design choice, and this
port does not reproduce it.

Because the C++ reference is not a usable oracle for circular mode, its
output is validated by properties instead of goldens: an independent
small-n enumeration (naive modular loops, no prefix sums or suffix-sum
walks) cross-checked against the production path; rotation invariance of
the partition function, every per-base signal, and the MFE location under a
rotation of the input sequence; bit-parity between every non-wrapped
window's energy in circular mode and the same window's energy in linear
mode; a closed-form check against an analytic homopolymer case; and
parallel equivalence/determinism across worker counts. The full validation
plan is in `docs/superpowers/specs/2026-11-03-circular-mode-design.md` §6.

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
port produces in linear mode — except perloop points computed with `--N
auto`, where this port deliberately diverges from the C++ (see Circular mode
above) — is checked against a patched build of the original C++ rlooper2 at
commit `455983bad1e9d49eb0e484f7746e60b16bd9868d`, at a tolerance of
`max(1e-9, 1e-5·|reference|)`. The goldens live under
`testdata/golden/` and can be regenerated at any time against a local
rlooper2 checkout with `scripts/gen-golden.sh`; the patches applied to the
C++ build (documented in `testdata/golden/README.md`) only remove undefined
behavior and do not change any computed value. Circular mode (`--circular`)
is validated differently, by the properties described in Circular mode
above, because the C++ reference's own circular implementation is broken
and is not a usable oracle.

A handful of rlooper2 behaviors that affect output values are preserved
deliberately rather than "fixed": the peak-calling threshold default above,
and signals that never credit the last base of each candidate window (an
off-by-one between a window's energy length and the length used for
residual/signal bookkeeping in the original). The nucleation (junction) free
energy `--a` is likewise folded into the ground-state energy only
(`ground state = K·α²/2 − a`) and excluded from every other structure's
reported energy — so the per-structure energies in `--dump` and `--top`
output are raw supercoiling-plus-base-pairing energies, not values relative
to the ground state. There are four places this port departs from rlooper2
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
port does not reproduce. Fourth, `--N auto`'s effectiveness under `perloop`:
the C++ `simulation_B` (perloop's equivalent) silently ignored
`auto_domain_size`, untested behavior rather than a deliberate choice, while
this port makes `--N auto` (and the circular `--N` default) effective under
`perloop` in both circular and linear runs (see Circular mode above).

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

`--circular` roughly doubles the window count for the same sequence (every
start has `n-minlength` lengths available instead of up to `n-start`), so
runtime and memory scale accordingly: on `res/pfc53_full.fa` (3,906nt),
`--dump` peaks at roughly 1.0GB RSS and writes a 344MB dump file under
`--circular`, versus roughly 0.5GB RSS and a 183MB dump file in linear mode.
