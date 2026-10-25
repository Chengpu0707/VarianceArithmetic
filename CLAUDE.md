# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Shape

Variance arithmetic: a floating-point type carrying a `(value, variance)` pair, described in the paper http://arxiv.org/abs/2410.01223. The same algorithms are implemented three times — C++ (`Cpp/`), Java (`Java/`), Python (`Python/`) — and the implementations must produce identical results.

The LaTeX paper (`Latex/VarianceArithmetic_Short.tex`, `Latex/VarianceArithmetic_Deduct.tex`) is the source of truth: code is checked against the tex, and `README.md` maps each paper section to its code, tests, and analysis notebooks. When changing an algorithm, change all three languages and keep them consistent with the paper.

Data flow: unit tests in each language write `*.txt` analysis files to `{Cpp,Java,Python}/Output/` → Jupyter notebooks in `IPyNb/` plot them → figures are saved as `*.png`/`*.pdf` under `Latex/` for the paper. `Maxima/*.wxmx` holds symbolic derivations for the matrix section.

**Before the first run**, create `Cpp/Output/`, `Cpp/bin/`, `Java/Output/`, and `Python/Output/` (all gitignored).

## Build and Test

### C++ (`Cpp/`)

Header-only library; each `Test*.cpp` is a standalone executable with its own `main()`. No project file. g++ is MSYS2 MinGW (`C:\msys64\mingw64\bin\g++.exe`, set in `.vscode/tasks.json`).

```
# single test (from Cpp/); VS Code's build task uses the same output dir plus -g -O0 -lstdc++exp
g++ -std=c++23 -o bin/TestXxx.exe TestXxx.cpp && ./bin/TestXxx.exe

# every Test*.cpp with -O3 -march=native, pass/fail summary (note: writes *.exe into Cpp/, not bin/)
./run_all_tests.sh
```

All code must compile and pass under `-std=c++03`, `c++11`, `c++17`, `c++20`, and `c++23`:
```
for std in c++03 c++11 c++17 c++20 c++23; do
  g++ -std=$std -o bin/TestXxx.exe TestXxx.cpp && ./bin/TestXxx.exe
done
```

On Windows, run executables from PowerShell (`& ".\bin\TestXxx.exe"`) to avoid MSYS2 bash segfaults.

### Java (`Java/`)

No Maven/Gradle. Sources and JUnit 4 tests live together in `Java/src/{Type,Func,Stats}` (packages `Type`, `Func`, `Stats`); jars are in `Java/lib`. VS Code's Java extension builds to `Java/bin` (see `.vscode/settings.json`). From the command line, **run from the repo root** (tests write to `./Java/Output/`):

```
CP="Java/lib/commons-math3-3.6.1.jar;Java/lib/junit-4.13.2.jar;Java/lib/hamcrest-core-1.3.jar"
javac -d Java/bin -cp "$CP" $(find Java/src -name '*.java')
java -cp "Java/bin;$CP" org.junit.runner.JUnitCore Type.TestVarDbl
```

### Python (`Python/`)

Interpreter is `.venv/Scripts/python.exe` (needs `numpy`, `scipy`, `sympy`). `pytest.ini` points pytest at `Python/` with `test*.py` files. Run from the repo root:

```
.venv/Scripts/python.exe -m pytest                                      # all
.venv/Scripts/python.exe -m pytest Python/testMoment.py::TestXxx::test_yyy   # one test
```

### Gated, slow, and ordered tests

- `Python/testManual.py` and `Python/testCompare.py` have `SKIP_TEST = True` at the top. Flip to `False` only to run them locally; **always commit as `True`**.
- `testManual.TestNormal` / `TestUniform` generate the samples that `testManual.TestBoundingRange` reads, so run them first.
- `testCompare.py` diffs selected C++/Java/Python `Output/` files and must run **last**, after all three languages' tests have populated their outputs.
- FFT full sweeps are very slow (C++ ~18 h unoptimized / ~3 h with `-O3 -march=native` via `.\bin\TestFFT.exe Test`; Java `TestFFT.dump_Order_2_19()` ~2 h; Python `testManual.Test_FFT_Order` ~1 week). They **resume** from an existing `Output/FFT_2_19.txt` — delete it for a clean run.
- The C++ and Java matrix drivers likewise skip `(size, noise)` cells already present in an existing `AdjMatrix_*.txt`.
- Tests use random sampling and occasionally fail (e.g. C++ `TestStat`); rerun before investigating.

## C++ Architecture

Header include graph (arrows point to what a header includes):

```
ulp.h ◄── Test.h ◄── VarDbl.h ◄─┬── Taylor.h ──► Moment.h
                                ├── Interval.h
Stat.h ◄────────────────────────┼── IndexSin.h
                                └── Matrix.h
FFT.h ──► IndexSin.h, Interval.h, Taylor.h, Stat.h
```

Note that the core type `VarDbl.h` depends on the test framework `Test.h`, because `VarDbl` has a built-in `assertEqual()` that calls `test::assertAlmostEqual`.

- `ulp.h`: units-in-the-last-place computation, with C++03-compatible type traits (`_is_floating_point_c03`, `_is_integral_c03`).
- `VarDbl.h`: the `VarDbl` value/variance type with operator overloads.
- `Moment.h`: bounded moments of Normal/Uniform distributions used by the Taylor expansion.
- `Taylor.h`: statistical Taylor expansion for `sin`, `exp`, `log`, `pow`, and polynomials. Its exception hierarchy (`TaylorIdException` + derived) needs `virtual ~XxxException() throw() {}` for C++03, because `std::runtime_error` has a `throw()` destructor.
- `Stat.h`: `Stat<T>` and `Histogram<T>` counters.
- `Interval.h`: interval type used by FFT.
- `IndexSin.h` / `FFT.h`: precomputed sine tables and the variance-arithmetic FFT.
- `Matrix.h`: Laplace-expansion determinant/adjugate; `runMatrixAnalysis()` is driven by `RunMatrixAnalysis.cpp` (`BenchMatrix.cpp` benchmarks it).
- `Test.h`, `TestTaylor.h`, `TestVarDbl.h`, `TestFFT.h`: the home-grown test framework and shared test helpers.

C++03-through-C++23 compatibility patterns used throughout:
- `VDBL_TMPL`, `STAT_TMPL`, `HISTO_TMPL` macros hide C++20 `requires` clauses.
- `#if __cplusplus >= 201103L` (and `202302L` etc.) guards around lambdas, `nullptr`, initializer lists, `<random>`, `<chrono>`, `<stacktrace>`.
- `> >` spacing in nested templates.
- `_C03FloatTag` / `_C03IntTag` tag structs for constructor SFINAE disambiguation.
- `_Optional<T>` stands in for `std::optional` below C++17.
- Unscoped enums accessed through the enclosing class (`IndexSin::Quart`, not `IndexSin::SinSource::Quart`).

## Java and Python Notes

- Java `VarDbl` has no operator overloading; arithmetic uses methods, with `*InPlace()` variants.
- Python `analytic.py` provides solutions for arbitrary analytic functions; Python also has extras with no C++/Java port (`movingLineFit.py`, `recursiveSin.py`, `hydrogenSpectra.py`).
