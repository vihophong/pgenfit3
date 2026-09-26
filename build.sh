#!/usr/bin/env bash
# Build the pgenfit3 model + fitter (and optionally the simulators) for one
# decay network -- or for a MIXTURE of several implanted species -- into
# BUILD_DIR. fit_decay_curve must be run from BUILD_DIR (it reads the path
# files there).
#
# Usage: ./build.sh <parmsex_file>[,<parmsex_file2>,...] <build_dir> [withsim]
#   one file : path.txt + decayModel_cal.cc (builder buildDecayModel)
#   several  : species i is generated in sp<i>/ with prefix s<i>_ (builder
#              buildDecayModel_s<i>_, path file path_s<i>_.txt); all are
#              linked into one libdecayModel.so
#   withsim  : also build simulation and simulation_mix
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <parmsex_file>[,<parmsex_file2>,...] <build_dir> [withsim]" >&2
    exit 1
fi
IFS=',' read -r -a PARMS_IN <<< "$1"
PARMS=()
for p in "${PARMS_IN[@]}"; do
    [ -f "$p" ] || { echo "parmsex file not found: $p" >&2; exit 1; }
    PARMS+=("$(cd "$(dirname "$p")" && pwd)/$(basename "$p")")
done
BUILD_DIR="$2"
WITHSIM="${3:-}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

CFLAGS=$(root-config --cflags)
LIBS=$(root-config --libs)
RPATH="-Wl,-rpath,$(root-config --libdir)"
CXX="g++ -std=c++17 -O2 -fsized-deallocation"

$CXX $CFLAGS -I"$SCRIPT_DIR" "$SCRIPT_DIR/decaypath.cc" "$SCRIPT_DIR/main.cc" -o main $LIBS $RPATH

MODEL_SRCS=()
if [ ${#PARMS[@]} -eq 1 ]; then
    ./main "${PARMS[0]}" > main.log
    MODEL_SRCS+=(decayModel_cal.cc)
    cat > decayModelRegistry.cc <<'EOF'
#include <decayModel.hh>
std::vector<DecayModelBuilder> decayModelBuilders(){ return {buildDecayModel}; }
EOF
else
    DECLS=""; LIST=""
    for i in "${!PARMS[@]}"; do
        pfx="s${i}_"
        mkdir -p "sp$i"
        (cd "sp$i" && ../main "${PARMS[$i]}" "$pfx" > main.log)
        cp "sp$i/path.txt" "path_${pfx}.txt"
        MODEL_SRCS+=("sp$i/decayModel_cal.cc")
        DECLS+="DecayModel* buildDecayModel_${pfx}(RooArgList&, RooRealVar&, RooRealVar&);"$'\n'
        LIST+="${LIST:+,}buildDecayModel_${pfx}"
    done
    printf '#include <decayModel.hh>\n%s' "$DECLS" > decayModelRegistry.cc
    echo "std::vector<DecayModelBuilder> decayModelBuilders(){ return {$LIST}; }" >> decayModelRegistry.cc
fi
$CXX -shared -fPIC $CFLAGS -I"$SCRIPT_DIR" "${MODEL_SRCS[@]}" decayModelRegistry.cc -o libdecayModel.so \
    $LIBS -lRooFit -lRooFitCore $RPATH
$CXX $CFLAGS -I"$SCRIPT_DIR" "$SCRIPT_DIR/fit_decay_curve.cxx" -o fit_decay_curve \
    -L. -ldecayModel $LIBS -lRooFit -lRooFitCore -Wl,-rpath,'$ORIGIN' $RPATH
if [ "$WITHSIM" = "withsim" ]; then
    $CXX $CFLAGS -I"$SCRIPT_DIR" "$SCRIPT_DIR/decaypath.cc" "$SCRIPT_DIR/mainsimulation.cc" \
        "$SCRIPT_DIR/simulation.cc" -o simulation $LIBS $RPATH
    $CXX $CFLAGS -I"$SCRIPT_DIR" "$SCRIPT_DIR/decaypath.cc" "$SCRIPT_DIR/mainsimulation_mix.cc" \
        "$SCRIPT_DIR/simulation.cc" -o simulation_mix $LIBS $RPATH
fi
echo "built in $BUILD_DIR (${#PARMS[@]} implanted species)"
