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

By default a state carries one to four simultaneous denatured runs, which is SIST's own convention and its published structural cap. `--max-runs N` (1 to 4) narrows that cap. `--max-runs 1` restricts the ensemble to the fully base-paired ground state plus every single melted window — the model this tool computed before multi-run support was added — and its output is unchanged from that earlier version, byte for byte.

In multi-run mode the length limits depend on how many runs a state has. Every run of every state is between `--minlength` and `--maxlength` bases long (defaults 1 and 250). On top of that, the two runs of a two-run state may total no more than `--maxlength` bases between them, and the runs of a three- or four-run state may total no more than 200 bases between them. That 200 is a constant in SIST's source and `--maxlength` does not move it, so raising `--maxlength` widens only the one- and two-run states, while lowering it narrows all of them. Both totals are reached at default settings on ordinary sequences.

```
# SIST-format report on stdout (multi-run, up to four runs per state)
golooper-sidd -f plasmid.fa --circular

# Single-run model only
golooper-sidd -f plasmid.fa --circular --max-runs 1

# Write the report to a file and also write wig tracks
golooper-sidd -f sequence.fa --clamp -o report.txt --tracks results
```

`--tracks results` writes `results_sidd_bpprob.wig` (P(x)) and `results_sidd_gx.wig` (G(x)). In multi-run mode the report is SIST's full line set: a parameter and ensemble header, the one-, two-, three- and four-run state counts, the total number of states, the partition function, the scaling factor, the average number of runs, the transition probability and the melted-base count, followed by a `Position  Base  P(x)  G(x)` table (SIST's run time line is not printed). The report does not record `--max-runs`, so runs that differ only in a non-binding cap produce identical output, and a zero three- or four-run count means none were found at these parameters, not that those classes were not modelled. With `--max-runs 1` the report is a shorter layout, matching what SIST prints when only its single-run search is active.

SIST's recommended input range is 1,500-10,000 bp. Multi-run mode is far slower than single-run because its state space is far larger, as in SIST itself: measured on one machine, about 3 s at 4.4 kb, about 70 s at 10 kb and a few minutes at 20 kb, against well under a second for `--max-runs 1` on the same inputs (memory roughly 60 MB at 4.4 kb to 270 MB at 20 kb). A long multi-run job has not hung. Much longer sequences at the default sigma make every single-domain state's Boltzmann weight underflow, and the tool reports an error rather than printing wrong numbers. `--sigma` must be less than or equal to 0; SIST's positive-magnitude `-s` value corresponds to a negative `--sigma` here.

### Relationship to SIST

The calculations follow [SIST](https://github.com/CCPBioSim/SIST) (MIT license), read at commit `cc206d0`.

Two enumeration conventions here differ deliberately from SIST's, and both apply only to states with more than one melted run.

- **Repeated run lengths.** When a state contains two or more runs of the same length, this implementation counts each unordered set of runs once. SIST's multi-run searches visit such a set under both orderings.
- **The retention window.** A state is kept if its energy is within `--threshold` of the lowest energy in the whole ensemble, which is found in a first pass before any state is accumulated. SIST advances that window as its search proceeds, so the set of retained states depends on the order in which they are visited.

These are the only intended departures, and they do change the numbers. On pBR322 with `--clamp` at default settings, compared with an unmodified SIST build: state counts are roughly 1% to 4% lower, the total number of states 2.0% lower, the partition function 1.6% lower, the melted-base count 0.04% higher, and the largest per-base probability difference is 4.2% relative. Constants, energetics, geometry and report format mirror SIST throughout.

**What is validated, and against what.**

- **Single-run mode** (`--max-runs 1`) is byte-identical to a build of SIST whose two-, three- and four-run searches are switched off, over eight recorded outputs in `testdata/golden/sidd/`. `scripts/gen-sidd-goldens.sh` regenerates them and records the build details.
- **Multi-run mode** is compared numerically rather than byte for byte, over nine recorded outputs in `testdata/golden/sidd-multirun/`, regenerated by `scripts/gen-multirun-goldens.sh`. Counts and labels must agree exactly, and printed floating-point values to a relative 1e-6. The build used for that comparison adopts the same repeated-run-length convention described above, so the two programs enumerate the same states and only their arithmetic can differ. Those sums run over millions of states, so the multi-run accumulators use Neumaier compensated summation to bound the rounding error that plain left-to-right addition accumulates; single-run mode keeps plain summation, which is what allows its byte-for-byte match. The retention-window convention is left as SIST has it; instead each case is admitted only after an instrumented build confirms that no state encountered later in the search lowers the minimum energy, which makes the two window conventions coincide for that case.
- **Four-run states** are not covered by any of these comparisons. A parameter sweep found no setting that both produces them and satisfies that admission check within a workable runtime, and the sweep was not exhaustive. They are checked instead against an independent brute-force enumerator, which agrees with the engine to a relative 1e-12 on short sequences.

`testdata/golden/sidd/README.md` and `testdata/golden/sidd-multirun/README.md` record how each set of reference outputs was produced, including the SIST commit, the build, and the checks each case had to pass.

SIST at these settings sums over states with up to four melted runs, and at pBR322 defaults its ensemble averages about 1.94 runs per state. The single-run model is therefore a physical approximation of it: on pBR322 at clamp defaults, SIST reports 79.5644 melted bases against 43.3706 for `--max-runs 1`, and the weaker site near position 3250 reaches a probability of 0.915 in SIST and 0.109 in the single-run model. Competition between R-loops and bubbles is not implemented.

If you use `golooper-sidd`, cite the SIST paper, which its authors require for any use of these algorithms, and the paper describing the strand-separation method:

> Zhabinskaya, D., Madden, S., & Benham, C. J. (2015). SIST: stress-induced structural transitions in superhelical DNA. *Bioinformatics*, 31(3), 421-422.
>
> Fye, R. M. and Benham, C. J. (1999). Exact method for numerically analyzing a model of local denaturation in superhelically stressed DNA. *Physical Review E*, 59, 3408-3426.

(Citation strings are from the SIST README, `CCPBioSim/SIST` at `cc206d0`.)

SIST's MIT license and the list of material derived from SIST are in [NOTICE](NOTICE).

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
| `--max-runs` | none | maximum simultaneous denatured runs per state, 1 to 4, default 4; SIST has no such parameter |
| `--N` | none | superhelical domain size; see below |
| `--minlength`, `--maxlength` | none | bounds on the length of a single melted run, defaults 1 and 250. In multi-run mode `--maxlength` also caps the combined length of a two-run state, while three- and four-run states are capped at a combined 200 bp that it does not change; SIST has no flag for either bound (its maximum window size is fixed at 250) |
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
