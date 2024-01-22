#include <cassert>
#include <cmath> 
#include <iostream>
#include <limits>
#include <stacktrace>   // not supported in gcc 13.2.0
#include <string>


#ifndef __Test_h__
#define __Test_h__

struct Test {
    // unit in the last place
    static double ulp(double x);

    // print current stack trace and exit
    static void fail(std::string msg = "");
    static void assertTrue(bool expression, std::string msg = "");
    static void assertEquals(double x, double y, double delta = 0, std::string msg = "");
        // ulp comparison when delta == 0
    template<typename T> static void assertEqual(const T& x, const T& y, std::string msg = "");
};


inline double Test::ulp(double x) 
{
    if (x > 0)
        return std::nexttoward(x, std::numeric_limits<double>::infinity()) - x;
    else 
        return x - std::nexttoward(x, -std::numeric_limits<double>::infinity());
}

inline void Test::fail(std::string msg) {
    if (!msg.empty())
        std::cout << msg << '\n';
#ifdef __cpp_lib_stacktrace
     std::cout << std::stacktrace::current() << '\n';
#endif
    assertTrue(false);
}


inline void Test::assertTrue(bool expression, std::string msg) {
    if (expression) 
        return;
    if (!msg.empty())
        std::cout << msg << '\n';
#ifdef __cpp_lib_stacktrace
    std::cout << std::stacktrace::current() << '\n';
#endif
    assert(expression);
}


inline void Test::assertEquals(double x, double y, double delta, std::string msg)
{
    if (delta == 0)
        delta = ulp(x);
    if (((x - delta) <= y) && (y <= (x + delta)))
        return;
    if (!msg.empty())
        std::cout << msg << '\n';
    std::cout << std::scientific << y << " != " << x << "+-" << delta << std::endl;
#ifdef __cpp_lib_stacktrace
    std::cout << std::stacktrace::current() << std::endl;
#endif
    assert(((x - delta) <= y) && (y <= (x + delta)));
}

template<typename T> 
inline void Test::assertEqual(const T& x, const T& y, std::string msg)
{
    if (x == y)
        return;
    if (!msg.empty())
        std::cout << msg << '\n';
    std::cout << std::scientific << y << " != " << x << std::endl;
#ifdef __cpp_lib_stacktrace
    std::cout << std::stacktrace::current() << std::endl;
#endif
    assert(x == y);
}

#endif  // __Test_h__