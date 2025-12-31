# VarianceArithmetic

## Purpose

A floating-point arithmetic with value and uncertainty (variance) pair.
The arithmetic is described in Latex/VarianceArithmetic,pdf.
It is also published at http://arxiv.org/abs/2410.01223.
The code and analysis associated with the paper will be pointed out here.

## Target

The code is targeted for Vs Code IDE.
For the ease of sharing the source code, no project file is provided.
Instead .vscode/settings.json and .vscode/tasks.json needs to be modified for each client's enviornment, such as the location of g++:
```
{
    "tasks": [
        {
            "type": "cppbuild",
            "label": "C/C++: g++.exe build active file",
            "command": "C:\\msys64\\mingw64\\bin\\g++.exe",
```


## Implementing Languages and Folders

The C++ implementation is under ./Cpp, as *.h files, for easier souce code sharing.
The cpp files contains tests using a self-made simple unit test framework provided by test.h.
There is no project file so each cpp file should be complied to be executed by itself.
Run all unit test will generate txt files under ./Cpp/Output

The Java implementation is under ./Java, containing both source code and test code.
The needed 3 jar libraries is under ./Java/lib.
Run all unit test will generate txt files under ./Java/Output

The Python implementation is under ./Python, containing both source code and test code.
Run all unit test will generate txt files under ./Python/Output
Specifically, analytic.py allows the solution to any analytic functions.
Not included is the python interpreter and libraries (numpy, scipy, sympy) which needs to be downloaded and set up, such as by the python pip tool. 

The txt files under ./{Cpp, Java, Python}/Output are analyzed by Jupyter Notebook code ipynb under ./IPyNb/as *.ipynb files.
Different implementation should give the same result.

The figures in *.ipynb files are saved as *.png files under ./Latex, and converted to *.pdf files under ./Latex before they are used in the ./Latex/VarianceArithmetic.tex file.
The Latex style file is ./Latex/intmacros.sty


# Unit Tests

All unit tests should pass.
Before running unit tests, please make sub folder ./{Cpp, Java, Python}/Output.

## Manual test

For Python, the manual tests in Python/testManual.py needs to run first, because Python/testStatBounding.py depends on TestNormal and TestUniform in Python/testManual.py.
To enable manual test in Python:
 1) Turn SKIP_TEST to False.  It should be committed as True.
 2) Run manually all the test labelled  @unittest.skipIf(SKIP_TEST, '...') from the command prompt as:
``` 
    ...\VarianceArithmetic\Python> Python -m unittest testManual.TestNormal
```

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
For other compilation options, please change .vscode/tasks.json and rerun the C++ FFT test.
The full FFT test took about 3 hours with the following compilation flag: -O3 march=native

For Python, run the manual command:
``` 
    ...\VarianceArithmetic\Python> Python -m unittest testManual.Test_FFT_Order.test_Gaussian
    ...\VarianceArithmetic\Python> Python -m unittest testManual.Test_FFT_Order.test_White
```
The full FFT test took about a week.


## Last Unit tests to run

testCompare.py compares the results of C++, Java, and Python for the selected tests, and it should be run last.

IPyNb/ExeTime.ipynb comoares the executioin speed of C++, Java, and Python.


## Outlier

Outliers is applied only when in FFT calculation 
 * the sine source is library, and
 * the signal is either sine or cosine without added noise, and 
 * the test is reverse transformation, and 
 * the expected value is 0+/-0 and
 * the value is error deviation whose average value is 1, and
 * An outlier minimal threshold of 1e14 is applied, to filter out 1e-16+/-1e-32 vs 0~0, such as:
```
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=32769, normaliized error outlier 8604408562923685.0 between 6.123234e-17~7.116e-33 and 0.000000e+00~0.000e+00
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=65536, normaliized error outlier -8604408562923685.0 between 0.000000e+00~0.000e+00 and 1.224647e-16~1.423e-32
For signal=Sin freq=1 noiseType=Gaussian noise=0 test=Reverse index=98305, normaliized error outlier -8604408562923685.0 between -6.123234e-17~7.116e-33 and 0.000000e+00~0.000e+00
```

# Failure of unit test

Because of using random variable in the test, occassionally, some test can fail due to sampling.
Usually rerun the test will pass.


# Fuction Description and Related files

The Jupiter files in ./IpyNb folder shows that difference implmentation of variance arithmetic has identical results.

./Python/testCompare.py also compares selected files from the Output folders from difference implmentation directly to show no difference in the results.


## Statistical Tools

A statistical counter, and a hostogram counter only for Deviation == 1:
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

Calculation of bounded momentum for Normal and Uniform distributions in _Sub Section Bounding Asymptote_ of _Section Statistical Taylor Expansion_:
 * Cpp/momentum.h
 * Java/Types/Momentum.java: calculate uncertainty as well.
 * Python/Momentum.py

Analysis:
 * IPyNb/Momentum.ipynb: with linear fit.

The unit tests:
 * Cpp/TestMomentum.cpp
 * Java/Types/TestMomentum.java
 * Python/testMomentum.py


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
 * Java/Types/TestVarDbl.java, Java/Types/TestVarDblAdd.java, Java/Types/TestVarDblMutiply.java
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
 * Python/testPolynominal.cpp

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

Anaylysis:
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


## _Section: Regressive Generation of Sin and Cos_

Code:
 * Python/regressiveSin.py

Unit tests:
 * Python/testRegressiveSin.py

Analysis:
 * IPyNb/RegressiveSin.ipynb


# Unit Test Coverage

## Python

Please download: 
```
$ pip install coverage
'''

Run all the test*.py files together, which will take much longer (~3930 secons) than the normal unit test (~820 seconds).
'''
VarianceArithmetic\Python> coverage run -m unittest testCompare.py testFFT.py testHisto.py testIndexSin.py testMatrix.py testMomentum.py testMovingLineFit.py testPolynominal.py testRegressiveSin.py testSinCos.py testStatBounding.py testTaylor.py testValDbl.py
Ran 256 tests in 3930.139s
OK (skipped=8)
```

Display the result:
```
VarianceArithmetic\Python> coverage report -m                                          
Name                   Stmts   Miss  Cover   Missing                                                                               
----------------------------------------------------
fft.py                   692    100    86%   163-165, 215, 217, 255-256, 259-262, 269-272, 315-316, 328-329, 355-356, 358-359, 368, 392-394, 399, 403, 412, 416, 420, 423, 444, 450-452, 457, 461, 464, 472, 501, 523-564, 587, 597, 605, 613, 617, 621, 632, 682, 686, 695, 699, 761-762, 776-777, 799, 818-819, 826-827, 864-866
histo.py                  86      6    93%   20, 33, 37, 53, 56, 102
indexSin.py              198     17    91%   42, 80-81, 125, 133-134, 158-161, 166, 170, 173, 177, 180, 183, 212, 232
matrix.py                125     11    91%   70, 80, 82, 90, 95, 106, 109, 159-161, 182
momentum.py              184     84    54%   17, 22, 27, 31, 43-57, 61-77, 88, 97, 100, 102, 104, 111, 115-157, 173, 178
movingLineFit.py          45      2    96%   68-69
regressiveSin.py          75      9    88%   38, 40, 42, 46, 58, 77, 84, 92, 98
taylor.py                306     50    84%   21, 66, 82-102, 156, 159, 163, 243-247, 263-265, 270-271, 277-279, 288-289, 293-295, 331-333, 350-351, 357, 377
varDbl.py                192     10    95%   59, 75-76, 167-168, 200, 248, 260-261, 275
'''

## C++

In .vscode/tasks.json, each file should be compiled with the following options:
```
"-ftest-coverage",
"-fprofile-arcs",
```
These options will slow down the execution.

'''
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

File 'C:/Users/Cheng/OneDrive/Documents/Github/VarianceArithmetic/Cpp/Momentum.h'
Lines executed:100.00% of 31
Creating 'Momentum.h.gcov'

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
'''