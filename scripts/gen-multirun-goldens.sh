#!/usr/bin/env bash
set -euo pipefail

FORCE=0
case "${1:-}" in
"") ;;
--force) FORCE=1 ;;
*)
    echo "usage: $0 [--force]" >&2
    exit 2
    ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SIST_DIR="$(cd "${SIST_DIR:-$REPO_ROOT/../SIST}" && pwd)"
DEDUP_PATCH="$SCRIPT_DIR/sist-dedup-fix.patch"
PROBE_PATCH="$SCRIPT_DIR/sist-probe.patch"
FIXTURES="$REPO_ROOT/testdata/fixtures"
GOLDEN_DIR="$REPO_ROOT/testdata/golden/sidd-multirun"
EXPECTED_COMMIT="cc206d0d22b0c631e8ee801e8cbd94aa61adbf90"

SIST_COMMIT="$(git -C "$SIST_DIR" rev-parse HEAD)"
if [ "$SIST_COMMIT" != "$EXPECTED_COMMIT" ]; then
    echo "SIST is at $SIST_COMMIT, expected $EXPECTED_COMMIT" >&2
    exit 1
fi

REFERENCE="$SIST_DIR/tests/regression/reference/v1.0.0/melting.txt"
TOY_FA="$FIXTURES/pbr322.toy.fa"
SYNTH_FA="$FIXTURES/sidd_synth1600.fa"
ROT_FA="$FIXTURES/pbr322_rot310.toy.fa"
for f in "$DEDUP_PATCH" "$PROBE_PATCH" "$REFERENCE" "$TOY_FA" "$SYNTH_FA" "$ROT_FA"; do
    if [ ! -f "$f" ]; then
        echo "missing file: $f" >&2
        exit 1
    fi
done
if ! cmp -s "$TOY_FA" "$SIST_DIR/tests/regression/data/pbr322.toy.fa"; then
    echo "$TOY_FA differs from the SIST copy" >&2
    exit 1
fi

BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT
STAGE_DIR="$BUILD_DIR/staged"
mkdir -p "$STAGE_DIR"

python3 -I -c "
s=''.join(l.strip() for l in open('$TOY_FA'))
r=s[310:]+s[:310]
[print(r[i:i+70]) for i in range(0,len(r),70)]" > "$BUILD_DIR/rot310.expected.fa"
if ! cmp -s "$ROT_FA" "$BUILD_DIR/rot310.expected.fa"; then
    echo "$ROT_FA is not the 310 bp rotation of $TOY_FA" >&2
    exit 1
fi

STOCK_DIR="$BUILD_DIR/stock"
DEDUP_DIR="$BUILD_DIR/dedup"
PROBE_DIR="$BUILD_DIR/probe"
cp -R "$SIST_DIR/src/trans_three" "$STOCK_DIR"
cp -R "$SIST_DIR/src/trans_three" "$DEDUP_DIR"
cp -R "$SIST_DIR/src/trans_three" "$PROBE_DIR"

echo "--- building unpatched qsidd ---"
make -C "$STOCK_DIR"

echo "--- build sanity: unpatched melting run vs reference ---"
SANITY_DIR="$BUILD_DIR/sanity"
mkdir "$SANITY_DIR"
cp "$TOY_FA" "$SANITY_DIR/pbr322.toy.fa"
(cd "$SANITY_DIR" && "$STOCK_DIR/qsidd" -b -p -r -f pbr322.toy.fa > stdout.txt)
if ! diff <(grep -v '^Run time' "$SANITY_DIR/stdout.txt") <(grep -v '^Run time' "$REFERENCE") > /dev/null; then
    echo "unpatched build does not reproduce $REFERENCE" >&2
    exit 1
fi

echo "--- applying sist-dedup-fix.patch (patch -p1) and building DEDUP ---"
(cd "$DEDUP_DIR" && patch -p1 < "$DEDUP_PATCH")
make -C "$DEDUP_DIR"

echo "--- applying sist-dedup-fix.patch and sist-probe.patch and building PROBE ---"
(cd "$PROBE_DIR" && patch -p1 < "$DEDUP_PATCH" && patch -p1 < "$PROBE_PATCH")
make -C "$PROBE_DIR"

run_case() {
    local case_name="$1"
    local fasta="$2"
    shift 2
    local run_dir="$BUILD_DIR/run/$case_name"
    local fasta_name start end elapsed improvements
    fasta_name="$(basename "$fasta")"
    mkdir -p "$run_dir"
    cp "$fasta" "$run_dir/$fasta_name"
    echo "--- running case: $case_name ---"
    start="$(date +%s)"
    (cd "$run_dir" && "$DEDUP_DIR/qsidd" -b -p -r "$@" -f "$fasta_name" > stdout.raw)
    end="$(date +%s)"
    elapsed=$((end - start))
    (cd "$run_dir" && "$PROBE_DIR/qsidd" -b -p -r "$@" -f "$fasta_name" > probe.raw 2> probe.err)
    if ! diff <(grep -v '^Run time' "$run_dir/stdout.raw") <(grep -v '^Run time' "$run_dir/probe.raw") > /dev/null; then
        echo "$case_name: PROBE stdout differs from DEDUP stdout" >&2
        exit 1
    fi
    improvements="$(grep -c 'old_min_E=' "$run_dir/probe.err" || true)"
    if [ "$improvements" != "0" ]; then
        echo "$case_name: PROBE gate failed, min_E improved after the 1R pre-pass:" >&2
        grep 'old_min_E=' "$run_dir/probe.err" >&2
        exit 1
    fi
    if [ "$(grep -c '^PROBE level=.R start ' "$run_dir/probe.err")" != "3" ]; then
        echo "$case_name: PROBE did not report the three level-start lines" >&2
        exit 1
    fi
    mkdir -p "$STAGE_DIR/$case_name"
    grep -v '^Run time' "$run_dir/stdout.raw" > "$STAGE_DIR/$case_name/stdout.txt"
    echo "$case_name: probe gate passed, DEDUP run time ${elapsed}s"
}

verify_case() {
    local case_name="$1"
    local expected_rows="$2"
    local kind="$3"
    local golden="$STAGE_DIR/$case_name/stdout.txt"
    local rows
    if [ ! -s "$golden" ]; then
        echo "missing or empty: $golden" >&2
        exit 1
    fi
    for level in one two three four; do
        if [ "$(grep -Ec "^Number of $level-run states = " "$golden")" != "1" ]; then
            echo "$golden lacks exactly one 'Number of $level-run states' line" >&2
            exit 1
        fi
    done
    rows="$(grep -Ec $'^[0-9]+\t[ACGT]\t' "$golden")"
    if [ "$rows" != "$expected_rows" ]; then
        echo "$golden has $rows table rows, expected $expected_rows" >&2
        exit 1
    fi
    case "$kind" in
    degenerate)
        for level in two three four; do
            if ! grep -Eq "^Number of $level-run states = 0$" "$golden"; then
                echo "$golden is marked degenerate but has nonzero $level-run states" >&2
                exit 1
            fi
        done
        ;;
    populated | hot)
        if grep -Eq '^Number of two-run states = 0$' "$golden"; then
            echo "$golden is marked $kind but has zero two-run states" >&2
            exit 1
        fi
        if [ "$kind" = hot ] && grep -Eq '^Number of three-run states = 0$' "$golden"; then
            echo "$golden is marked hot but has zero three-run states" >&2
            exit 1
        fi
        ;;
    *)
        echo "unknown case kind: $kind" >&2
        exit 1
        ;;
    esac
}

run_case mr-clamp-copoly-defaults "$TOY_FA"
run_case mr-clamp-nn-defaults "$TOY_FA" -n
run_case mr-circular-copoly-defaults "$TOY_FA" -c
run_case mr-circular-nn-synth "$SYNTH_FA" -c -n
run_case mr-clamp-copoly-alt "$TOY_FA" -T 315 -i 0.1 -s 0.04 -t 10
run_case mr-clamp-copoly-hot "$TOY_FA" -T 320 -t 15 -s 0.04
run_case mr-circular-nn-defaults "$TOY_FA" -c -n
run_case mr-clamp-copoly-synth "$SYNTH_FA"
run_case mr-origin-straddle "$ROT_FA" -c

echo "--- verifying goldens ---"
verify_case mr-clamp-copoly-defaults 4361 populated
verify_case mr-clamp-nn-defaults 4361 populated
verify_case mr-circular-copoly-defaults 4361 populated
verify_case mr-circular-nn-synth 1600 populated
verify_case mr-clamp-copoly-alt 4361 degenerate
verify_case mr-clamp-copoly-hot 4361 hot
verify_case mr-circular-nn-defaults 4361 populated
verify_case mr-clamp-copoly-synth 1600 populated
verify_case mr-origin-straddle 4361 populated

echo "--- checking staged goldens against existing goldens ---"
changed=0
for staged in "$STAGE_DIR"/*/stdout.txt; do
    case_name="$(basename "$(dirname "$staged")")"
    existing="$GOLDEN_DIR/$case_name/stdout.txt"
    if [ -f "$existing" ] && ! cmp -s "$staged" "$existing"; then
        echo "$case_name: staged golden differs from the committed golden:" >&2
        diff "$existing" "$staged" | head -20 >&2 || true
        changed=1
    fi
done
if [ "$changed" = "1" ] && [ "$FORCE" != "1" ]; then
    echo "refusing to overwrite differing goldens; rerun with --force to replace them" >&2
    exit 1
fi

echo "--- installing goldens ---"
mkdir -p "$GOLDEN_DIR"
for staged in "$STAGE_DIR"/*/stdout.txt; do
    case_name="$(basename "$(dirname "$staged")")"
    mkdir -p "$GOLDEN_DIR/$case_name"
    cp "$staged" "$GOLDEN_DIR/$case_name/stdout.txt"
done

echo "--- writing README ---"
DEDUP_SHA="$(shasum -a 256 "$DEDUP_PATCH" | cut -d' ' -f1)"
PROBE_SHA="$(shasum -a 256 "$PROBE_PATCH" | cut -d' ' -f1)"
COMPILER="$(c++ --version | head -1)"
{
    cat <<README
# SIST multi-run golden reference data

Generated by \`scripts/gen-multirun-goldens.sh\`:

\`\`\`
SIST_DIR=/path/to/SIST bash scripts/gen-multirun-goldens.sh
\`\`\`

- SIST commit: \`$SIST_COMMIT\`
- Compiler: \`$COMPILER\`
- \`scripts/sist-dedup-fix.patch\` sha256: \`$DEDUP_SHA\`
- \`scripts/sist-probe.patch\` sha256: \`$PROBE_SHA\`
- \`SIST_DIR\` environment variable: path to the SIST checkout, default \`../SIST\`
  relative to the repo root. That default is wrong from a nested worktree
  (\`.claude/worktrees/<name>\`); set \`SIST_DIR\` explicitly there.

## Build

\`src/trans_three\` is copied three times to scratch directories and built with plain
\`make\` (\`c++ -O2 -Wall\`): stock, stock plus the dedup patch (DEDUP), and stock plus
the dedup and probe patches (PROBE).

## Oracle

In SIST's 2R/3R/4R searches a pair of equal-length runs is visited under both
orderings (\`n1 == n2\` restarts the inner iterator at \`begin()\` plus one rather
than after the outer iterator), so such a pair enters the sums twice.
\`scripts/sist-dedup-fix.patch\` (\`patch -p1\` inside a scratch copy of
\`src/trans_three\`) starts the inner iterator after the outer one, matching the
convention golooper-sidd uses, so that the two programs enumerate the same
states. The goldens are the output of that build, with \`^Run time\` lines removed.

Before patching, the stock binary is run as \`./qsidd -b -p -r -f pbr322.toy.fa\`
and must match \`tests/regression/reference/v1.0.0/melting.txt\` (ignoring
\`^Run time\`); the script aborts otherwise.

## Soundness gate (order-dependent theta window)

SIST lowers \`min_E\` (and the \`max_E = min_E + theta\` window) as the 2R/3R/4R
searches discover lower-energy states, so the set of enumerated states depends on
search order. A golden is valid only if no \`min_E\` improvement occurs after the
1R pre-pass. \`scripts/sist-probe.patch\` (applied on top of the dedup patch in a
second scratch copy) writes to stderr one \`PROBE level=NR start\` line at the start of
each level and one \`old_min_E=\` line on any improvement. Per case, the script
requires: PROBE stdout identical to DEDUP stdout (ignoring \`^Run time\`), zero
\`old_min_E=\` lines, and exactly three level-start lines. A case failing the gate
aborts the script.

## Cases

All cases run \`qsidd -b -p -r <flags> -f <fasta>\` in a directory holding a copy of
the fixture; \`-f\` must be last.

| case | fixture | flags | multi-run expectation |
|---|---|---|---|
README
    for entry in \
        "mr-clamp-copoly-defaults|pbr322.toy.fa|(none)|populated" \
        "mr-clamp-nn-defaults|pbr322.toy.fa|-n|populated" \
        "mr-circular-copoly-defaults|pbr322.toy.fa|-c|populated" \
        "mr-circular-nn-synth|sidd_synth1600.fa|-c -n|populated" \
        "mr-clamp-copoly-alt|pbr322.toy.fa|-T 315 -i 0.1 -s 0.04 -t 10|degenerate (all multi-run counts zero)" \
        "mr-clamp-copoly-hot|pbr322.toy.fa|-T 320 -t 15 -s 0.04|populated, three-run states > 0" \
        "mr-circular-nn-defaults|pbr322.toy.fa|-c -n|populated" \
        "mr-clamp-copoly-synth|sidd_synth1600.fa|(none)|populated" \
        "mr-origin-straddle|pbr322_rot310.toy.fa|-c|populated"; do
        IFS='|' read -r name fasta flags expectation <<< "$entry"
        echo "| \`$name\` | \`$fasta\` | \`$flags\` | $expectation |"
    done
    cat <<README

Every listed case passed the probe gate. \`mr-clamp-copoly-alt\` is intentionally
degenerate: at these parameters SIST finds no two-, three- or four-run states, so it
pins the multi-run line set printing zero counts. \`verify_case\` enforces each
expectation.

## v1 parameter-set coverage

The v1 matrix (\`testdata/golden/sidd/README.md\`) has eight parameter sets. This
matrix covers seven: clamp-copoly-defaults, clamp-nn-defaults, circular-copoly-defaults,
circular-nn-defaults, clamp-copoly-synth, circular-nn-synth and clamp-copoly-alt
(\`mr-clamp-copoly-alt\`). \`circular-nn-alt\` (sidd_synth1600.fa, \`-c -n -T 300 -i 0.05
-s 0.07 -t 11\`) is not covered. \`mr-clamp-copoly-hot\` has no v1 counterpart. The
\`mr-circular-nn-defaults\` and \`mr-clamp-copoly-synth\` cases close audit finding 6a
for the two sets the first matrix dropped.

## Origin-straddle case

\`mr-origin-straddle\` runs the circular copolymeric defaults on \`pbr322_rot310.toy.fa\`,
pbr322 rotated by 310 bp so the dominant melting site (1-based positions 286-329 in
the original golden tables, as printed) straddles the origin: it spans the
last and first bases of the sequence, so its runs cross index 0. The M13 separation
rule is not rotation-invariant (the three-versus-one gap asymmetry depends on runs
touching the last index), so its counts may legitimately differ from \`mr-circular-copoly-defaults\`; the
golden comparison is the gate, and no rotation-invariance is asserted.

## No four-run case

There is no \`mr-clamp-copoly-4r\` case: no point found by the 2026-11-04 parameter
sweep (theta, stress, temperature and salt inside SIST's no-WARNING envelope) gave
\`Number of four-run states\` > 0 with the probe gate passing within the 10-minute
run-time eligibility bound. The decisive points on \`pbr322.toy.fa\`:

| flags | four-run states | outcome |
|---|---|---|
| \`-t 15 -s 0.07\` | 41134 | 24 post-1R min_E drops, gate fail |
| \`-t 15 -s 0.07 -i 0.0001\` | 1882 | gate fail, at 71 s |
| \`-t 15 -s 0.05 -T 320\` | 573 | 17 min, over the bound, gate not run |

Unexplored regions: stress 0.05 with T between 310 and 320, and salt between 0.02 and
0.05.

## Fixtures

- \`testdata/fixtures/pbr322.toy.fa\`: verbatim copy of \`SIST/tests/regression/data/pbr322.toy.fa\`.
- \`testdata/fixtures/sidd_synth1600.fa\`: see \`testdata/golden/sidd/README.md\`.
- \`testdata/fixtures/pbr322_rot310.toy.fa\`: pbr322.toy.fa rotated by 310 bp, generated
  from the repo root with:

\`\`\`
python3 -I -c "
s=''.join(l.strip() for l in open('testdata/fixtures/pbr322.toy.fa'))
r=s[310:]+s[:310]
[print(r[i:i+70]) for i in range(0,len(r),70)]" > testdata/fixtures/pbr322_rot310.toy.fa
\`\`\`
README
} > "$GOLDEN_DIR/README.md"

echo "--- writing MANIFEST.sha256 ---"
(
    cd "$REPO_ROOT"
    shasum -a 256 scripts/sist-dedup-fix.patch scripts/sist-probe.patch testdata/golden/sidd-multirun/*/stdout.txt
) > "$GOLDEN_DIR/MANIFEST.sha256"

echo "done. goldens written to $GOLDEN_DIR"
