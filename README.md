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
 * Python/matrix.py

Unit test:
 * Python/testMatrix.py
 * Python/testManual.py: class TestAdjugateManually

Analysis:
 * IPyNb/AdjugateMatrix.ipynb
 * IPyNb/MatrixCondition.ipynb
 * Maxima/Matrix_2.wxmx, Maxima/Matrix_2.wxmx


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