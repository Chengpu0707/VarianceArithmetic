/*
A minimal unittest tool box
*/

#include <cassert>
#include <iostream>
#include <sstream>
#include <stacktrace>   // not supported in gcc 13.2.0
#include <string>

#include "ulp.h"

#ifndef __Test_h__
#define __Test_h__
namespace test 
{

struct AssertException : public std::runtime_error 
{
    explicit AssertException(const std::string& what_arg) : runtime_error(what_arg) 
    {
#ifdef __cpp_lib_stacktrace
        std::cout << std::stacktrace::current() << '\n';
#endif
        std::cout << what() << '\n';
    }
};


// throw AssertException, which will print current stack trace and exit if not caught
static void fail(std::string msg = "");
static void assertTrue(bool expression, std::string msg = "");
static void assertFalse(bool expression, std::string msg = "");
static void assertEquals(double x, double y, double delta = 0, std::string msg = "");
    // ulp comparison when delta == 0
template<typename T, typename U> static void assertEqual(const T& x, const U& y, std::string msg = "");


inline void fail(std::string msg) {
    if (!msg.empty())
        std::cout << msg << '\n';
    std::ostringstream os;
    os << "fail(" << msg << ")";
    throw AssertException(os.str());
}


inline void assertTrue(bool expression, std::string msg) 
{
    if (expression) 
        return;
    std::ostringstream os;
    os << "assertTrue(" << msg << ")";
    throw AssertException(os.str());
}

inline void assertFalse(bool expression, std::string msg) 
{
    return assertTrue(!expression, msg);
}


inline void assertEquals(double x, double y, double delta, std::string msg)
{
    if (std::isfinite(x) != std::isfinite(y)) {
        std::ostringstream os;
        os << "assertEquals("<< std::scientific << y << " != " << x;
        if (!msg.empty())
            os << ", " << msg;
        os << ")";
        throw AssertException(os.str());
    }
    if (!std::isfinite(x))
        return;
    if (delta == 0)
        delta = var_dbl::ulp(x);
    if (((x - delta) <= y) && (y <= (x + delta)))
        return;
    std::ostringstream os;
    os << "assertEquals("<< std::scientific << y << " != " << x << "+-" << delta;
    if (!msg.empty())
        os << ", " << msg;
    os << ")";
    throw AssertException(os.str());
}

template<typename T, typename U> 
inline void assertEqual(const T& x, const U& y, std::string msg)
{
    if (x == y)
        return;
    std::ostringstream os;
    os << "assertEqual(" << std::scientific << y << " != " << x;
    if (!msg.empty())
        os << ", " << msg;
    os << ")";
    throw AssertException(os.str());
}

} // namespace test
#endif  // __Test_h__