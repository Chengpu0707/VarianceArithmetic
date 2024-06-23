/*
A minimal unittest tool box
*/

#include <cmath>
#include <cassert>
#include <iostream>
#include <ranges>
#include <sstream>
#include <stacktrace>
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
        std::cout << what() << '\n';
#if _HAS_CXX23
        std::stacktrace st;
        std::cout << st << '\n';
#endif
    }
};


// throw AssertException, which will print current stack trace and exit if not caught
static void fail(std::string msg = "");
static void assertTrue(bool expression, std::string msg = "");
static void assertFalse(bool expression, std::string msg = "");
static void assertAlmostEqual(long double x, long double y, long double delta = 0, std::string msg = "");
    // ulp comparison when delta == 0
template<typename T, typename U> static void assertEqual(const T& x, const U& y, std::string msg = "");
    // generic value comparison
template<typename T, typename U> static void assertEquals(const T& sX, const U& sY, std::string msg = "");
    // generic collection comparison


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


inline void assertAlmostEqual(long double x, long double y, long double delta, std::string msg)
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
    if (((y - delta) <= x) && (x <= (y + delta)))
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

template<typename T, typename U> 
inline void assertEquals(const T& sX, const U& sY, std::string msg)
{
    std::ostringstream os;
    if (sX.size() != sY.size()) {
        os << "assertEquals(size: " << sY.size() << " != " << sX.size();
        if (!msg.empty())
            os << ", " << msg;
        os << ")";
        throw AssertException(os.str());
    }
    size_t i = 0;
    for (auto [x, y]: std::ranges::views::zip(sX, sY)) {
        if (x != y) {
            os << "assertEquals(" << i << ": " << sX[i] << " != " << sY[i];
            if (!msg.empty())
                os << ", " << msg;
            os << ")";
            throw AssertException(os.str());
        }
        ++i;
    }
}

} // namespace test
#endif  // __Test_h__