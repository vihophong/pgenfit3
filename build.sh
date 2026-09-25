#!/usr/bin/env bash
# Build the pgenfit3 model + fitter (and optionally the simulator) for one
# decay network into BUILD_DIR. Everything a fit needs (path.txt,
# libdecompModel.so, fit_decay_curve) ends up there; fit_decay_curve must
# be run from BUILD_DIR since it reads ./path.txt.
#
# Usage: ./build.sh <parmsex_file> <build_dir> [withsim]
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <parmsex_file> <build_dir> [withsim]" >&2
    exit 1
fi
PARMSEX="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
BUILD_DIR="$2"
WITHSIM="${3:-}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

CFLAGS=$(root-config --cflags)
LIBS=$(root-config --libs)
RPATH="-Wl,-rpath,$(root-config --libdir)"
CXX="g++ -std=c++17 -O2 -fsized-deallocation"

$CXX $CFLAGS -I"$SCRIPT_DIR" "$SCRIPT_DIR/decaypath.cc" "$SCRIPT_DIR/main.cc" -o main $LIBS $RPATH
./main "$PARMSEX" > main.log
$CXX -shared -fPIC $CFLAGS -I"$SCRIPT_DIR" decayModel_cal.cc -o libdecayModel.so \
    $LIBS -lRooFit -lRooFitCore $RPATH
$CXX $CFLAGS -I"$SCRIPT_DIR" "$SCRIPT_DIR/fit_decay_curve.cxx" -o fit_decay_curve \
    -L. -ldecayModel $LIBS -lRooFit -lRooFitCore -Wl,-rpath,'$ORIGIN' $RPATH
if [ "$WITHSIM" = "withsim" ]; then
    $CXX $CFLAGS -I"$SCRIPT_DIR" "$SCRIPT_DIR/decaypath.cc" "$SCRIPT_DIR/mainsimulation.cc" \
        "$SCRIPT_DIR/simulation.cc" -o simulation $LIBS $RPATH
fi
echo "built in $BUILD_DIR ($(head -1 path.txt) species)"
