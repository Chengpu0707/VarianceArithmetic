# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test

All C++ tests are compiled and run from the `Cpp/` directory. Each `.cpp` file is a standalone test executable — there is no project file.

**Build a single test (from `Cpp/`):**
```
g++ -std=c++23 -o TestXxx.exe TestXxx.cpp && ./TestXxx.exe
```

**Build across all supported standards:**
```
for std in c++03 c++11 c++17 c++20 c++23; do
  g++ -std=$std -o TestXxx.exe TestXxx.cpp && ./TestXxx.exe
done
```

On Windows, use PowerShell and `& ".\TestXxx.exe"` to run executables (avoids bash segfault issues with MSYS2).

**Required before first run:** create `Cpp/Output/` directory.

**FFT test** is extremely slow unoptimized (~18 hours). Use `-O3 -march=native` (~3 hours). Full test: `.\TestFFT.exe Test`.

**Flaky tests:** `TestStat` (testWhite/testGaussian) uses live random seeds and can occasionally fail due to sampling. Rerun to pass — this is expected and documented.

The g++ path is configured in `.vscode/tasks.json` (`C:\msys64\mingw64\bin\g++.exe`).

## Architecture

**Language:** C++ header-only library. All library code is in `*.h` files under `Cpp/`. Each `TestXxx.cpp` includes the headers it needs and has its own `main()`.

**Core dependency chain:**
```
ulp.h → VarDbl.h → Momentum.h → Taylor.h → IndexSin.h / FFT.h
```

- `ulp.h`: ULP (units in the last place) computation with C++03-compatible type traits (`_is_floating_point_c03`, `_is_integral_c03`, `_C03FloatTag`/`_C03IntTag` tag structs for SFINAE disambiguation)
- `VarDbl.h`: Core `VarDbl` type — a `(value, variance)` pair with operator overloads. Uses `VDBL_TMPL` macro to collapse C++20 `requires` clauses vs plain `template<typename T>`
- `Momentum.h`: Bounded moments for Normal/Uniform distributions used in Taylor expansion
- `Taylor.h`: Statistical Taylor expansion for `sin`, `exp`, `log`, `pow`, and polynomial functions. Exception hierarchy (`TaylorIdException` and 5 derived) must use `virtual ~XxxException() throw() {}` in C++03 due to `std::runtime_error`'s `throw()` dtor
- `Stat.h`: `Stat<T>` and `Histogram<T>` counters. Uses `STAT_TMPL`/`HISTO_TMPL` macros; `_Optional<T>` polyfill for C++17 `std::optional`
- `IndexSin.h`: Precomputed sine/cosine tables for FFT precision
- `FFT.h`: Variance-arithmetic FFT implementation
- `Test.h` / `TestTaylor.h`: Lightweight test framework (`validate_func`, `search_edge`, `stat_func`)

**C++ standard compatibility:** All files target C++03 through C++23. Key patterns used:
- `STAT_TMPL`, `HISTO_TMPL`, `VDBL_TMPL` macros abstract away C++20 `requires` clauses
- `#if __cplusplus >= 201103L` guards for lambdas, `nullptr`, initializer lists, `<random>`, `<chrono>`, etc.
- Template `> >` spacing (not `>>`) for C++03 nested templates
- `_C03FloatTag`/`_C03IntTag` tag structs for constructor SFINAE disambiguation in C++03
- `_Optional<T>` struct replaces `std::optional` below C++17
- Unscoped enum access via enclosing class (`IndexSin::Quart`) not `IndexSin::SinSource::Quart`

**Output:** Tests write analysis files to `Cpp/Output/` (e.g., `SinEdge.txt`, `ExpStat.txt`). These are consumed by Jupyter notebooks in `IPyNb/` for visualization.

**Other implementations:** Java (`Java/`) and Python (`Python/`) implement the same algorithms. `Python/testCompare.py` validates that all three produce identical output files.
