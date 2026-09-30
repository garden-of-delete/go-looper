#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RLOOPER2_DIR="$(cd "$REPO_ROOT/../rlooper2" && pwd)"

FIXTURE="$REPO_ROOT/testdata/fixtures/synthetic200.fa"
PFC53="$REPO_ROOT/res/pfc53_full.fa"
GOLDEN_DIR="$REPO_ROOT/testdata/golden"

if [ ! -f "$FIXTURE" ]; then
    echo "missing fixture: $FIXTURE" >&2
    exit 1
fi
if [ ! -f "$PFC53" ]; then
    echo "missing fixture: $PFC53" >&2
    exit 1
fi

if ! git -C "$RLOOPER2_DIR" diff --quiet HEAD --; then
    echo "$RLOOPER2_DIR has uncommitted changes to tracked files; refusing to build from a dirty HEAD" >&2
    exit 1
fi
RLOOPER2_COMMIT="$(git -C "$RLOOPER2_DIR" rev-parse HEAD)"

BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT

SRC_DIR="$BUILD_DIR/rlooper2"
mkdir -p "$SRC_DIR"
git -C "$RLOOPER2_DIR" archive "$RLOOPER2_COMMIT" | tar -x -C "$SRC_DIR"

EQ_MODEL="$SRC_DIR/Rloop_equilibrium_model.cpp"
EQ_MODEL_LINE192="$(sed -n '192p' "$EQ_MODEL")"
if [ "$EQ_MODEL_LINE192" != '    std::vector<char>::iterator b_0,b_1;' ]; then
    echo "Rloop_equilibrium_model.cpp:192 does not match the expected text; refusing to patch" >&2
    echo "found: $EQ_MODEL_LINE192" >&2
    exit 1
fi
perl -i -pe 's/^    std::vector<char>::iterator b_0,b_1;$/    std::vector<char>::iterator b_0,b_1=start;/ if $. == 192' "$EQ_MODEL"

SIM_CPP="$SRC_DIR/simulation.cpp"
SEED_LINE_COUNT="$(grep -c '^    seed = 0;$' "$SIM_CPP")"
if [ "$SEED_LINE_COUNT" != "1" ]; then
    echo "expected exactly one 'seed = 0;' line in simulation.cpp, found $SEED_LINE_COUNT" >&2
    exit 1
fi
perl -i -pe 's/^(    seed = 0;)$/$1\n    bedfile = false;\n    residuals = false;/' "$SIM_CPP"

MAKEFILE="$SRC_DIR/Makefile"
MAKEFILE_LINE4="$(sed -n '4p' "$MAKEFILE")"
if [ "$MAKEFILE_LINE4" != 'CFLAGS = -g -std=c++11' ]; then
    echo "Makefile:4 does not match the expected text; refusing to patch" >&2
    echo "found: $MAKEFILE_LINE4" >&2
    exit 1
fi
perl -i -pe 's/^CFLAGS = -g -std=c\+\+11$/CFLAGS = -g -std=c++11 -O2/ if $. == 4' "$MAKEFILE"

echo "--- verifying patches landed ---"
grep -n 'b_1=start' "$EQ_MODEL"
grep -n 'bedfile = false;\|residuals = false;' "$SIM_CPP"
grep -n '^CFLAGS' "$MAKEFILE"

OBJS="gene.o model.o Rloop_equilibrium_model.o simulation.o structure.o windower.o Rloop_dynamic_model.o"

echo "--- compiling shared objects ---"
mkdir -p "$SRC_DIR/bin"
make -C "$SRC_DIR" $OBJS ensemble_analyzer_main.o

echo "--- linking rlooper ---"
(cd "$SRC_DIR" && g++ -g -std=c++11 -O2 -o bin/rlooper $OBJS ensemble_analyzer_main.o)

echo "--- building perloop ---"
(cd "$SRC_DIR" && g++ -g -std=c++11 -O2 -c perloop_main.cpp -o perloop_main.o)
(cd "$SRC_DIR" && g++ -g -std=c++11 -O2 -o bin/perloop $OBJS perloop_main.o)

RLOOPER_BIN="$SRC_DIR/bin/rlooper"
PERLOOP_BIN="$SRC_DIR/bin/perloop"
if [ ! -x "$RLOOPER_BIN" ] || [ ! -x "$PERLOOP_BIN" ]; then
    echo "expected binaries were not produced" >&2
    exit 1
fi

run_rlooper_case() {
    local case_name="$1"
    local infile="$2"
    shift 2
    local outdir="$SRC_DIR/golden/$case_name"
    mkdir -p "$outdir"
    echo "--- running case: $case_name ---"
    (cd "$outdir" && "$RLOOPER_BIN" "$infile" "out" "$@" > "stdout.txt")
}

run_rlooper_case defaults "$FIXTURE" --localaverageenergy --dump --residuals --top 5 --bedfile --sensitivity -1
run_rlooper_case sigma010 "$FIXTURE" --sigma -0.10
run_rlooper_case nauto "$FIXTURE" --N auto
run_rlooper_case homopolymer05 "$FIXTURE" --homopolymer 0.5
run_rlooper_case pfc53 "$PFC53"

echo "--- running case: perloop ---"
PERLOOP_OUTDIR="$SRC_DIR/golden/perloop"
mkdir -p "$PERLOOP_OUTDIR"
(cd "$PERLOOP_OUTDIR" && "$PERLOOP_BIN" "$FIXTURE" "out" -0.10 -0.02 > "stdout.txt")

echo "--- verifying outputs ---"
for case_name in defaults sigma010 nauto homopolymer05 pfc53; do
    outdir="$SRC_DIR/golden/$case_name"
    for suffix in _bpprob.wig _avgG.wig _mfe.wig _bpprob.bed _mfe.bed _extbpprob.wig; do
        f="$outdir/out$suffix"
        if [ ! -s "$f" ]; then
            echo "missing or empty: $f" >&2
            exit 1
        fi
    done
    stdout_file="$outdir/stdout.txt"
    if [ ! -s "$stdout_file" ]; then
        echo "missing or empty: $stdout_file" >&2
        exit 1
    fi
    if ! grep -q 'P(ground state)=' "$stdout_file"; then
        echo "$stdout_file does not contain 'P(ground state)='" >&2
        exit 1
    fi
    if grep -qi 'exception' "$stdout_file"; then
        echo "$stdout_file appears to contain exception text" >&2
        exit 1
    fi
done

DUMP_FILE="$SRC_DIR/golden/defaults/outSYNTH200_dump.txt"
if [ ! -s "$DUMP_FILE" ]; then
    echo "missing or empty: $DUMP_FILE" >&2
    exit 1
fi
DUMP_FIRST_ROW="$(sed -n '2p' "$DUMP_FILE")"
if ! echo "$DUMP_FIRST_ROW" | grep -Eq '^0 0 [^ ]+ 0$'; then
    echo "unexpected dump first data row: $DUMP_FIRST_ROW" >&2
    exit 1
fi

if ! grep -Eq '\.\.\.19900$' "$SRC_DIR/golden/defaults/stdout.txt"; then
    echo "expected a structure count of 19900 in defaults/stdout.txt" >&2
    exit 1
fi

if [ ! -s "$PERLOOP_OUTDIR/out" ]; then
    echo "missing or empty: $PERLOOP_OUTDIR/out" >&2
    exit 1
fi

echo "--- copying golden files into $GOLDEN_DIR ---"
for case_name in defaults sigma010 nauto homopolymer05 pfc53; do
    dest="$GOLDEN_DIR/$case_name"
    rm -rf "$dest"
    mkdir -p "$dest"
    cp "$SRC_DIR/golden/$case_name"/* "$dest/"
done

PERLOOP_DEST="$GOLDEN_DIR/perloop"
rm -rf "$PERLOOP_DEST"
mkdir -p "$PERLOOP_DEST"
cp "$PERLOOP_OUTDIR/out" "$PERLOOP_DEST/perloop_out.txt"
cp "$PERLOOP_OUTDIR/stdout.txt" "$PERLOOP_DEST/stdout.txt"

cat > "$GOLDEN_DIR/README.md" <<EOF
# C++ golden reference data

Generated by \`scripts/gen-golden.sh\` from a patched build of rlooper2.

rlooper2 commit: \`$RLOOPER2_COMMIT\`
Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)

## Patches applied to the build copy

1. \`Rloop_equilibrium_model.cpp:192\`: \`std::vector<char>::iterator b_0,b_1;\` -> \`std::vector<char>::iterator b_0,b_1=start;\` (fixes uninitialized-\`b_1\` UB in \`compute_structure\`; the other occurrence of this declaration, at line 218 in \`compute_external_structure\`, is untouched).
2. \`Simulation\` constructor (\`simulation.cpp\`): added \`bedfile = false;\` and \`residuals = false;\`.
3. \`Makefile:4\`: appended \`-O2\` to \`CFLAGS\`.

None of these change program behavior for the golden matrix below; they only remove undefined behavior from uninitialized fields/iterators.

## perloop output location

\`perloop_main.cpp\` takes \`sequence_file output_file superhelicity_lower_bound superhelicity_upper_bound\`. \`Simulation::simulation_B\` writes \`<superhelicity> <p_rloop>\` for each superhelicity step to the **file** passed as \`output_file\` (via \`ofstream outfile(argv[2], ios::out)\`); stdout only receives a progress print of the superhelicity value for each step, no \`P(ground state)=\` text. The golden reference is therefore the output **file**, copied to \`perloop/perloop_out.txt\`. \`perloop/stdout.txt\` is also kept, for reference only.

## Layout

- \`defaults/\`, \`sigma010/\`, \`nauto/\`, \`homopolymer05/\`, \`pfc53/\`: each holds the six \`simulation_A\` output files (\`out_bpprob.wig\`, \`out_avgG.wig\`, \`out_mfe.wig\`, \`out_bpprob.bed\`, \`out_mfe.bed\`, \`out_extbpprob.wig\`) plus \`stdout.txt\`. \`defaults/\` additionally has \`outSYNTH200_dump.txt\` (from \`--dump\`).
- \`perloop/\`: \`perloop_out.txt\` (the golden data) and \`stdout.txt\` (progress log only).
EOF

echo "done. golden files written to $GOLDEN_DIR"
