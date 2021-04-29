typedef long long bint;
#include <functional>
#include <algorithm>
#include <limits>
#include <ctime>
#include <cmath>
#include <string>
#include <thread>

// have NDEBUG undefined to evaluate assert(...)
//#define NDEBUG
#include <cassert>

// have ARMA_NO_DEBUG undefined to use boundary checks etc in armadillo
// #define ARMA_NO_DEBUG
#include <armadillo>

// undef possible stupid macros
#undef NAN
#undef min
#undef max
// end: undef possible stupid macros

#ifndef __EstEntSGconstants
#define __EstEntSGconstants
const double NaN = std::numeric_limits<double>::quiet_NaN();
const double Infinity = std::numeric_limits<double>::infinity();
#endif
