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

Right now, C++ and Java implementations contains less test cases such as for matrix calculations.
