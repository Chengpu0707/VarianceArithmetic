# VarianceArithmetic
A floating-point arithmetic with value and uncertainty (variance) pair.
The arithmetic is described in Latex/VarianceArithmetic
It is also published at http://arxiv.org/abs/2410.01223

The pyton implementation is under ./Python, containing both source code and test code.
Run all unit test will generate txt files under ./Python/Output
Specifically, analytic.py allows the solution to any analytic functions.

The C++ implementation is under ./Cpp, as *.h files, with a self-made simple unit test framework using assert.
Run all unit test will generate txt files under ./Cpp/Output

The Java implementation is under ./Java, containing both source code and test code.
Run all unit test will generate txt files under ./Java/Output

The txt files under ./{Python, Cpp, Java}/Output are analyzed by Jupyter Notebook code ipynb under ./IPyNb.
Different implementation should give the same result.

Some tests takes a long time to finish, so the test can resume from where leftoff in the result file.
For these tests, it is preferable to delete and remake the /Output folder, and manually run the test in a separated console.

Right now, C++ and Java implementations contains less test cases such as for matrix calculations.

Some performance numbrs in minutes on my old laptop:
Language  FFT/14 FFT/15 FFT/16 FFT/17  FFT/18 Adj/6 Adj/7 Adj/8
C++        17     36     76     141    300
Java
Python     92    272    463    1083   1756    11    118	  1425
