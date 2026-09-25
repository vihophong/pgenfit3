# Sourced by every pgenfit3 script. Needs ROOT >= 6.32 (RooExponential's
# negateCoefficient flag, RooFit::EvalBackend). If the ROOT on PATH is
# older, point PGENFIT3_ROOT at a newer install, e.g.
#   export PGENFIT3_ROOT=/data01/userdata/ROOT/root-6.40.02_patch
if [ -n "${PGENFIT3_ROOT:-}" ]; then
    export ROOTSYS="$PGENFIT3_ROOT"
    export PATH="$PGENFIT3_ROOT/bin:$PATH"
    export LD_LIBRARY_PATH="$PGENFIT3_ROOT/lib:${LD_LIBRARY_PATH:-}"
fi
if ! command -v root-config >/dev/null 2>&1; then
    echo "root-config not found; set PGENFIT3_ROOT to a ROOT >= 6.32 install" >&2
    exit 1
fi
_rv=$(root-config --version | tr '/' '.')
_rmaj=${_rv%%.*}; _rmin=${_rv#*.}; _rmin=${_rmin%%.*}
if [ "$_rmaj" -lt 6 ] || { [ "$_rmaj" -eq 6 ] && [ "$_rmin" -lt 32 ]; }; then
    echo "ROOT $_rv is too old (need >= 6.32); set PGENFIT3_ROOT to a newer install" >&2
    exit 1
fi
