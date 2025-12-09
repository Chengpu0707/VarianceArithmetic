# VarianceArithmetic

## Purpose

A floating-point arithmetic with value and uncertainty (variance) pair.
The arithmetic is described in Latex/VarianceArithmetic
It is also published at http://arxiv.org/abs/2410.01223

## Target

The code is targeted for Vs Code IDE.
For the ease of sharing source code, no project file is provided.

## Implementing Languages and Folders

The pyton implementation is under ./Python, containing both source code and test code.
Run all unit test will generate txt files under ./Python/Output
Specifically, analytic.py allows the solution to any analytic functions.
Not included is the python interpreter and libraries (numpy, scipy, sympy) which needs to be pip-ed. 

The C++ implementation is under ./Cpp, as *.h files, for easier souce code sharing.
The cpp files contains tests using a self-made simple unit test framework provided by test.h.
There is no project file so each cpp file should be complied to be executed by itself.
Run all unit test will generate txt files under ./Cpp/Output

The Java implementation is under ./Java, containing both source code and test code.
Run all unit test will generate txt files under ./Java/Output

The txt files under ./{Python, Cpp, Java}/Output are analyzed by Jupyter Notebook code ipynb under ./IPyNb.
Different implementation should give the same result.


# Unit Tests

## Manual test

The manual tests needs to run first before all the unit test can pass:
* Python: TestBoundingFactor relies on TestNormal and TestUniform
* Python: testCompare.py compares the results of Java, C++, and Python for the selected tests.

### FFT Test

The output file is ./{Java, Cpp, Python}/Output/FFT_2_19.txt.
If it exist, the test can be continue.
For a clean test, delete the file.

Java unit test TestFFT.dump_Order_2_19() takes about 2.1 hours

For debugging, C++ code is not compiled with any optimization, so it is much slower than java.
The full FFT test took 17.9 hours:
```
    ...\VarianceArithmetic\Cpp> .\TestFFT.exe Test
```

To enable manual test in Python:
 1) Turn SKIP_TEST to False.  It is committed as True.
 2) Run manually from the command prompt as:
``` 
    ...\VarianceArithmetic\Python> Python -m unittest testManual.Test_FFT_Order.test_Gaussian
```

## Outlier

Outliers is applied only when in FFT calculation 
 * the sine source is library, and
 * the signal is either sine or cosine without added noise, and 
 * the test is reverse transformation, and 
 * the expected value is 0~0 and
 * the value is error deviation whose average value is 1, and
 * An outlier minimal threshold of 1e14 is applied, to filter out 1e-16~1e-32 vs 0~0, such as:
```
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=32769, normaliized error outlier 8604408562923685.0 between 6.123234e-17~7.116e-33 and 0.000000e+00~0.000e+00
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=65536, normaliized error outlier -8604408562923685.0 between 0.000000e+00~0.000e+00 and 1.224647e-16~1.423e-32
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=98305, normaliized error outlier -8604408562923685.0 between -6.123234e-17~7.116e-33 and 0.000000e+00~0.000e+00
```


# Fuction Description and Related files

The Jupiter files in ./IpyNb folder shows that difference implmentation of variance arithmetic has identical results.

./Python/testCompare.py also compares selected files from the Output folders from difference implmentation directly to show no difference in the results.


## Bounded Momentum

Calculation of bounded momentum for Normal and Uniform distributions in _ Sub Section Bounding Asymptote _:
 * Cpp/momentum.h
 * Java/Types/Momentum.java: calculate uncertainty as well.
 * Python/Momentum.py
Analysis:
 * IPyNb/Momentum.ipynb: Linear fit.

The unit tests:
 * Cpp/TestMomentum.cpp
 * Java/Types/TestMomentum.java
 * Python/testMomentum.py

## VarDbl

The basic data type for variance arithmetic in _Sub Section Numerical Representation_ of  _Section Variance Arithmetic_:
 * Cpp/VarDbl.h
 * Java/Types/VarDbl.java: no operator override, with *InPlace() version for arithmetic operations.
 * Python/varDbl.py

The unit tests:
 * Cpp/TestVarDbl.cpp
 * Java/Types/TestVarDbl.java, Java/Types/TestVarDblAdd.java, Java/Types/TestVarDblMutiply.java
 * Python/testVarDbl.py

## Statistical Taylor Expansion

The numerical implementation of statistical Taylor expansion in _Section Variance Arithmetic_:
 * Cpp/Taylor.h
 * Java/Types/Taylor.java
 * Python/taylor.py

Analysis:
 * IPyNb/TaylorExpansion.ipynb: expansion of 1/(1 - x)

Direct application of statistical Taylor expansion in _Section Mathematical Library Functions_
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

Convergence edge calculation for x^c in _Sub Section Monotonic_ and _Sub Section Monotonic_ of _Section Variance Arithmetic_:
 * Cpp/TestExp.cpp, Cpp/TestLog.cpp, Cpp/TestSin.cpp, Cpp/TestPow.cpp
 * Java/TestConvergence.java
 * Python/testManual.py: class TestConvergence:

 Analysis:
  * IPyNb/ConvergeEdge.ipynb

## Polynomal Expansion

_Section Polynomial_
 * Cpp/TestPolynominal.cpp
 * Java/TestPolynominal.Java
 * Python/testPolynominal.cpp

Analysis: