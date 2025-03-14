#pragma once
#include <map>
#include <set>
#include "BSM.h"

// Normal CDF function (declaration)
double norm_cdf(double x);

// Normal PDF function
double norm_pdf(double x);

// Bisection Method
double bisection_method(BSM &bs, double market_price, bool debug = false);
