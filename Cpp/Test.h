#include <cassert>
#include <iostream>
#include <limits>
#include <stacktrace>
#include <string>


#ifndef __Test_h__
#define __Test_h__

struct Test {
    // unit in the last place
    static double ulp(double x);

    // print current stack trace and exit
    static void fail() {assertTrue(false); }
    static void assertTrue(bool expression);
    static void assertEquals(double x, double y, double delta = 0);
        // ulp comparison when  delta == 0
    template<typename T> static void assertEqual(const T& x, const T& y);
};


inline double Test::ulp(double x) 
{
    if (x > 0)
        return std::nexttoward(x, std::numeric_limits<double>::infinity()) - x;
    else 
        return x - std::nexttoward(x, -std::numeric_limits<double>::infinity());
}

inline void Test::assertTrue(bool expression) {
    if (expression) 
        return;
    //std::cout << std::stacktrace::current() << '\n';
    assert(expression);
}


inline void Test::assertEquals(double x, double y, double delta)
{
    if (delta == 0)
        delta = ulp(x);
    if (((x - delta) <= y) && (y <= (x + delta)))
        return;
    //std::cout << std::stacktrace::current() << '\n';
    assert(((x - delta) <= y) && (y <= (x + delta)));
}

template<typename T> inline void Test::assertEqual(const T& x, const T& y)
{
    if (x == y)
        return;
    //std::cout << std::stacktrace::current() << '\n';
    assert(x == y);
}

#endif  // __Test_h__