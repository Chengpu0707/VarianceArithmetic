#!/usr/bin/env bash
# Build and run every Test*.cpp in Cpp/. Each is a standalone executable
# (one main() per file). Output: per-test status + combined summary.
#
# Notes:
#   * TestFFT is slow (~3 hours under -O3 -march=native; ~18h unoptimized).
#     Run with no argument here for the short smoke test; pass `Test` for
#     the full sweep manually if needed.
#   * TestStat occasionally flakes on live RNG seeds.

set +e
cd "$(dirname "$0")"

CXX=g++
CXXFLAGS="-std=c++23 -O3 -march=native -I . -Wno-reorder"

pass=0
fail=0
sFail=()

for src in Test*.cpp; do
    name="${src%.cpp}"
    exe="${name}.exe"
    echo "=== building $src ==="
    if ! $CXX $CXXFLAGS -o "$exe" "$src" 2>&1 | tail -20; then
        echo "BUILD FAIL: $src"
        sFail+=("$name (build)")
        ((fail++))
        continue
    fi
    if [ ! -f "$exe" ]; then
        echo "BUILD FAIL (no exe produced): $src"
        sFail+=("$name (build)")
        ((fail++))
        continue
    fi
    echo "=== running $exe ==="
    start=$(date +%s)
    "./$exe" 2>&1 | tail -10
    rc=$?
    dt=$(( $(date +%s) - start ))
    if [ $rc -eq 0 ]; then
        echo "PASS: $name (${dt}s)"
        ((pass++))
    else
        echo "FAIL: $name (rc=$rc, ${dt}s)"
        sFail+=("$name (rc=$rc)")
        ((fail++))
    fi
done

echo
echo "=========================================="
echo "Summary: $pass passed, $fail failed"
if [ ${#sFail[@]} -gt 0 ]; then
    echo "Failures:"
    for f in "${sFail[@]}"; do
        echo "  - $f"
    done
fi
echo "=========================================="
exit $fail
