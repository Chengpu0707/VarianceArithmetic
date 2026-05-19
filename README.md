# VarianceArithmetic

## Purpose

A floating-point arithmetic with value and uncertainty (variance) pair.
The arithmetic is described as a PDF in http://arxiv.org/abs/2410.01223, and as a presentation in https://www.youtube.com/watch?v=dfSihfKt7Xc.
The code and analysis associated with the paper will be pointed out here.

## Target

The code is targeted for the VS Code IDE.
For ease of sharing the source code, no project file is provided.
Instead, `.vscode/settings.json` and `.vscode/tasks.json` need to be modified for each client's environment, such as the location of g++:
```
{
    "tasks": [
        {
            "type": "cppbuild",
            "label": "C/C++: g++.exe build active file",
            "command": "C:\\msys64\\mingw64\\bin\\g++.exe",
```


## Implementing Languages and Folders

The C++ implementation is under `./Cpp` as header (`*.h`) files for easier source-code sharing.
The `*.cpp` files contain tests using a self-made simple unit-test framework provided by `Test.h`.
There is no project file, so each `*.cpp` file should be compiled to be executed by itself.
Running all unit tests will generate `*.txt` files under `./Cpp/Output`.
All files are built and tested against C++23, C++20, C++17, C++11, and C++03 standards.

The Java implementation is under `./Java`, containing both source code and test code.
The 3 required jar libraries are under `./Java/lib`.
Running all unit tests will generate `*.txt` files under `./Java/Output`.

The Python implementation is under `./Python`, containing both source code and test code.
Running all unit tests will generate `*.txt` files under `./Python/Output`.
Specifically, `analytic.py` provides solutions for arbitrary analytic functions.
Not included are the Python interpreter and libraries (`numpy`, `scipy`, `sympy`), which need to be downloaded and set up (e.g., via `pip`).

The `*.txt` files under `./{Cpp, Java, Python}/Output` are analyzed by Jupyter notebooks under `./IPyNb/` as `*.ipynb` files.
Different implementations should produce the same result.

The figures in `*.ipynb` files are saved as `*.png` files under `./Latex`, and converted to `*.pdf` files under `./Latex` before being used in the `./Latex/VarianceArithmetic.tex` file.
The LaTeX style file is `./Latex/intmacros.sty`.

All codes are checked against the tex document using Claude.

# Unit Tests

All unit tests should pass.
Before running unit tests, please create the subfolder `./{Cpp, Java, Python}/Output`.

## Manual test

For Python, the manual tests in Python/testManual.py needs to run first, because Python/testStatBounding.py depends on TestNormal and TestUniform in Python/testManual.py.
To enable manual test in Python:
 1) Turn SKIP_TEST to False.  It should be committed as True.
 2) Run manually all the test labelled  @unittest.skipIf(SKIP_TEST, '...') from the command prompt as:
``` 
    ...\VarianceArithmetic\Python> Python -m unittest testManual.TestNormal
```

### FFT Test

The output file is `./{Java, Cpp, Python}/Output/FFT_2_19.txt`.
If it exists, the test will resume from where it left off.
For a clean test, delete the file.

The Java unit test `TestFFT.dump_Order_2_19()` takes about 2.1 hours.

For debugging, the C++ code is compiled without any optimization, so it is much slower than Java.
The full FFT test took 17.9 hours:
```
    ...\VarianceArithmetic\Cpp> .\TestFFT.exe Test
```
For other compilation options, change `.vscode/tasks.json` and rerun the C++ FFT test.
With `-O3 -march=native`, the full FFT test took about 3 hours.

For Python, run the manual command:
``` 
    ...\VarianceArithmetic\Python> Python -m unittest testManual.Test_FFT_Order.test_Gaussian
    ...\VarianceArithmetic\Python> Python -m unittest testManual.Test_FFT_Order.test_White
```
The full FFT test took about a week.


## Last Unit tests to run

`testCompare.py` compares the results of C++, Java, and Python for the selected tests, and should be run last.

`IPyNb/ExeTime.ipynb` compares the execution speed of C++, Java, and Python.


## Outlier

Outliers is applied only when in FFT calculation 
 * the sine source is library, and
 * the signal is either sine or cosine without added noise, and 
 * the test is reverse transformation, and 
 * the expected value is 0+/-0 and
 * the value is error deviation whose average value is 1, and
 * An outlier minimal threshold of 1e14 is applied, to filter out 1e-16+/-1e-32 vs 0~0, such as:
```
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=32769, normalized error outlier 8604408562923685.0 between 6.123234e-17~7.116e-33 and 0.000000e+00~0.000e+00
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=65536, normalized error outlier -8604408562923685.0 between 0.000000e+00~0.000e+00 and 1.224647e-16~1.423e-32
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=98305, normalized error outlier -8604408562923685.0 between -6.123234e-17~7.116e-33 and 0.000000e+00~0.000e+00
```

# Failure of unit test

Because the tests use random variables, occasionally a test can fail due to sampling.
Rerunning the test usually passes.


# Function Description and Related Files

The Jupyter notebooks in `./IPyNb/` show that the different implementations of variance arithmetic produce identical results.

`./Python/testCompare.py` also compares selected files from the `Output` folders of different implementations directly to show no difference in the results.


## Statistical Tools

A statistical counter, and a histogram counter only for Deviation == 1:
 * Cpp/Stats.h
 * Java/Stats/Stat.java, Java/Stats/Histogram.java
 * Python/histo.py

Unit test:
 * Cpp/TestStats.cpp
 * Java/Stats/TestStat.java, Java/Stats/TestHistogram.java
 * Python/testHisto.py


## _Section: Statistical Taylor Expansion_

### _Sub Section: Distributional Zero and Distributional Pole_

IpyNb/DistributionZeroAndPole.ipynb


### _Sub Section: Bounding Asymptote_

Calculation of bounded moment for Normal and Uniform distributions in _Sub Section Bounding Asymptote_ of _Section Statistical Taylor Expansion_:
 * Cpp/Moment.h
 * Java/Types/Moment.java: calculate uncertainty as well.
 * Python/moment.py

Analysis:
 * IPyNb/Moment.ipynb: with linear fit.

The unit tests:
 * Cpp/TestMoment.cpp
 * Java/Types/TestMoment.java
 * Python/testMoment.py


### _Sub Section: Statistical Bounding_ and  _Sub Section: Ideal Statistics_

Code as unit test:
 * TestNormal and TestUniform in Python/testManual.py
 * Python/testStatBounding.py

Analysis:
 * IPyNb/StatisticalBounding.ipynb


## _Section: Variance Arithmetic_

Code:
 * Cpp/Taylor.h
 * Java/Types/Taylor.java
 * Python/taylor.py

###  _Sub Section: Numerical Representation_

The basic data type for variance arithmetic in _Sub Section Numerical Representation_ of  _Section Variance Arithmetic_:
 * Cpp/VarDbl.h
 * Java/Types/VarDbl.java: no operator override, with *InPlace() version for arithmetic operations.
 * Python/varDbl.py

The unit tests:
 * Cpp/TestVarDbl.cpp
 * Java/Types/TestVarDbl.java, Java/Types/TestVarDblAdd.java, Java/Types/TestVarDblMultiply.java
 * Python/testVarDbl.py


### _Sub Section: Monotonic_ and _Sub Section: Positive_ 

Unit tests:
 * Cpp/TestExp.cpp, Cpp/TestLog.cpp, Cpp/TestSin.cpp, Cpp/TestPow.cpp
 * Java/TestConvergence.java
 * Python/testManual.py: class TestConvergence

 Analysis:
  * IPyNb/ConvergeEdge.ipynb


## _Section: Mathematical Library Functions_

Unit tests: 
 * Cpp/TestExp.cpp, Cpp/TestLog.cpp, Cpp/TestSin.cpp, Cpp/TestPow.cpp
 * Java/TestTaylor.java
 * Python/testTaylor.py

Analysis:
 * IPyNb/ExpStat{Cpp, Java, Python}.ipynb: for ./{Cpp, Java, Python}/Output/ExpStat.txt
 * IPyNb/LogStat{Cpp, Java, Python}.ipynb: for ./{Cpp, Java, Python}/Output/LogStat.txt
 * IPyNb/SinStat{Cpp, Java, Python}.ipynb: for ./{Cpp, Java, Python}/Output/SinStat.txt
 * IPyNb/PowStat{Cpp, Java, Python}.ipynb: for ./{Cpp, Java, Python}/Output/PowStat.txt
 * IPyNb/PowerError.ipynb: (x^c)^(1/c) - x error
 * IPyNb/ExpLogError.ipynb: log(e^x) - x and e^(log(x)) - x error.



## _Section: Polynomial_

Code:
 * Cpp/Taylor.h
 * Java/Types/Taylor.java
 * Python/taylor.py

Unit tests:
 * Cpp/TestPolynominal.cpp
 * Java/TestPolynominal.Java
 * Python/testPolynominal.py

Analysis:
 * IPyNb/TaylorExpansion.ipynb: expansion of 1/(1 - x)
 * IPyNb/Polynominal.ipynb
 * IPyNb/ConvergeEdge.ipynb


## _Section: Matrix Calculations_

Code:
 * Cpp/Matrix.h: bottom-up Laplace expansion (Formula 6.1/6.3) for `determ()` and `adjugate()`; `multiply()` for matrix product; static `runMatrixAnalysis()` that parallelizes the Hilbert + random-integer sweep via `std::async`.
 * Java/src/Func/Matrix.java: mirrors the C++ design (same bottom-up Laplace, same `multiply()` with inline-doubles fast path, same `runMatrixAnalysis(int minSize, int maxSize, int targetCells, double[] noises, String outDir, int threads)` static driver); `RandomMatrix` subclass populates noisy cells via `set()`.
 * Python/matrix.py

Unit test:
 * Cpp/TestMatrix.cpp: 15 standalone tests (construction, `get`/`set`, Laplace determinant, adjugate, the M·adj(M) = |M|·I identity, variance propagation, `RandomMatrix`).
 * Java/src/Func/TestMatrix.java: 16 JUnit tests — the same 15 plus `testAdjugateUncertaintyVsMonteCarlo` (a Monte-Carlo diagnostic that prints the Adj/Fwd/Rnd over- or under-prediction ratio at five noise levels, see explanation of the `Adj Norm Deviation ≈ 1 → drop` pattern below).
 * Python/testMatrix.py
 * Python/testManual.py: class TestAdjugateManually

Drivers and benchmarks:
 * Cpp/RunMatrixAnalysis.cpp: standalone executable wrapping `Matrix::runMatrixAnalysis`. CLI: `<minSize> <maxSize> <targetCells> <outDir>` (all positional, defaults `4 10 1000 Output`). Build: `g++ -std=c++23 -O2 -pthread -o RunMatrixAnalysis.exe RunMatrixAnalysis.cpp -lstdc++exp`.
 * Cpp/BenchMatrix.cpp: micro-benchmark for `calc()` at sizes 6/8/10; useful when regressing the inner-loop optimization (`std::map` → flat `vector<vector<VarDbl>>`, bit-trick set-bit iteration, in-place `+=`/`-=`, `__builtin_ctzll`/`popcountll`).
 * Java `Func.Matrix.main(String[] args)`: same CLI as the C++ driver, with optional 5th arg `threads` (default = `availableProcessors()`). Lower the thread count and pass `-Xmx...` to JVM if larger sizes OOM (sizes 14–15 need ≥ 8 GB / few threads on a 16 GB box).
 * Python/runMatrixAnalysis.py: original Python driver; the C++/Java versions are ports of it.

Output schema:
 * `MatrixCondition_{minSize}_{maxSize}.txt` — one row per matrix (Hilbert + Random samples). Columns: `Size, Type, Noise, Condition Number, Determinant Value, Determinant Uncertainty, Determinant Precision, Run Time`.
 * `AdjMatrix_{minSize}_{maxSize}.txt` — one row per `(size, noise)` cell, 58 columns: `Type, Noise, Size, Count` followed by 9 stat blocks for `{Adj, Fwd, Rnd} × {Unc, Val, Norm} × {Deviation, Mean, Minimum, Maximum, Count, Loss}`. `Adj` = noisy adjugate vs precise integer adjugate; `Fwd` = `M·adj − det·I`; `Rnd` = `M·adj − adj·M`.
 * The Java and C++ runners both **skip cells already present** in an existing `AdjMatrix_*.txt` (parsed at startup), so a previous run can be extended (e.g., rename `_4_10.txt` to `_4_16.txt` and rerun with the larger maxSize to fill in only sizes 10..15).

Analysis:
 * IPyNb/AdjugateMatrix.Python.ipynb: reads `Python/Output/AdjMatrix_4_9.txt`, plots 3D surfaces of the five error metrics (Adjugate Uncertainty/Error, Forward, Roundtrip, Multiple Error) plus normalized-error histograms and a log-noise power-law regression of Forward Error.
 * IPyNb/AdjugateMatrix.Cpp.ipynb: reads `Cpp/Output/AdjMatrix_4_16.txt`. Same surface plots adapted to the C++ schema (`Adj/Fwd/Rnd Unc/Val/Norm`), no histogram section since the C++/Java schema doesn't carry per-bucket counts. Includes the log-Forward power-law fit.
 * IPyNb/AdjMatrix.Java.ipynb: identical to `AdjugateMatrix.Cpp.ipynb` but reads `Java/Output/AdjMatrix_4_16.txt`.
 * IPyNb/MatrixCondition.ipynb
 * Maxima/Matrix_2.wxmx, Maxima/Matrix_2.wxmx

Notes on observed behaviors:
 * `Adj Norm Deviation` stays near 1 for noise ≲ 0.1 and drops rapidly (to ~0.5–0.6) at noise ≳ 1. The drop is the variance arithmetic over-predicting at large noise: Laplace's alternating-sign summation creates negative covariance between terms that `VarDbl.add` cannot see (it sums variances as if independent). `testAdjugateUncertaintyVsMonteCarlo` in `Java/src/Func/TestMatrix.java` confirms this with a Monte-Carlo comparison.
 * `Fwd Val Dev` at noise = 0, size ≥ 8 has heavy-tailed sampling variance — a few near-singular random matrices dominate the std. C++ and Java compute identical adj/M·adj values for any **single fixed** matrix; the orders-of-magnitude difference visible in the two AdjMatrix files reflects different random samples (each language has its own RNG, sample count is ~12 per cell at size 8). Fix: pass a deterministic seed or bump the sample count.


## _Section: Moving-Window Linear Regression_

Code:
 * movingLineFit.py

Unit test:
 * testMovingLineFit.py

Analysis:
 * IPyNb/MovingLineFit.ipynb


## _Section: FFT_

### _Sub Section: Modeling Errors of DFT_

IPyNb/FFT_unfaithful.ipynb


### _Sub Section: FFT_ and _Sub Section: Testing Signals_

Code:
 * Cpp/FFT.h, Cpp/TestFFT.h
 * Java/FFT.java
 * Python/fft.py

Unit tests:
 * Cpp/TestFFT.cpp
 * Java/TestFFT.java
 * Python/testFFT.py

Analysis:
 * IPyNb/FFT.{Cpp, Java, Python}.ipynb
 * IPyNb/FFT_Step_{Prec, Quart}.{Cpp, Java, Python}.ipynb
 * IPyNb/FFT_Step_Lib.SciPy.{Cpp, Java, Python}.ipynb
 * IPyNb/IndexSine_18.Quart.ipynb
 * IPyNb/SinDiff.ipynb
 

### _Sub Section: Trigonometric Library Errors_

Code:
 * Cpp/IndexSin.h
 * Java/Func/IndexSin.java
 * Python/indexSin.py

Unit tests:
 * Cpp/TestIndexSin.h
 * Java/Func/TestIndexSin.java
 * Python/testIndexSin.py

Analysis:
 * IPyNb/SinDiff.ipynb


## _Section: Recursive Generation of Sin and Cos_

Code:
 * Python/recursiveSin.py

Unit tests:
 * Python/testRecursiveSin.py

Analysis:
 * IPyNb/RecursiveSin.ipynb


# Unit Test Coverage

## Python and Java

Coverage for Python and Java is measured using the built-in coverage support in VS Code (via the Python and Java test extensions). Run the test suites from the VS Code Testing view with the "Run Tests with Coverage" action.

## C++

In .vscode/tasks.json, each file should be compiled with the following options:
```
"-ftest-coverage",
"-fprofile-arcs",
```
These options will slow down the execution.

```
...\VarianceArithmetic\Cpp> .\TestFFT.exe
...\VarianceArithmetic\Cpp> gcov .\TestFFT.cpp
PS C:\Users\Cheng\OneDrive\Documents\Github\VarianceArithmetic\Cpp> gcov .\TestFFT.cpp
File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/TestFFT.cpp'
Lines executed:86.26% of 444
Creating 'TestFFT.cpp.gcov'

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/IndexSin.h'
Lines executed:54.26% of 129
Creating 'IndexSin.h.gcov'

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/Taylor.h'
Lines executed:36.14% of 249
Creating 'Taylor.h.gcov'

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/Moment.h'
Lines executed:100.00% of 31
Creating 'Moment.h.gcov'

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/VarDbl.h'
Lines executed:84.87% of 119
Creating 'VarDbl.h.gcov'

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/VarDbl.h'
Lines executed:84.87% of 119
Creating 'VarDbl.h.gcov'

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/Stat.h'
Lines executed:100.00% of 92
Creating 'Stat.h.gcov'

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/Test.h'
Lines executed:29.23% of 65
Creating 'Test.h.gcov'

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/FFT.h'
Lines executed:100.00% of 74
Creating 'FFT.h.gcov'
```