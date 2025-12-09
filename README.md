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
Not included is the python interpreter and libraries (numpy, scipy, sympy) 
which needs to be piped by the donwloader. 

The C++ implementation is under ./Cpp, as *.h files, with a self-made simple unit test framework using assert.
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

