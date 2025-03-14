#include "util.h"
#include <cmath>
#include <numbers>
#include <array>

// Implementation of the Normal CDF function
double norm_cdf(double x)
{
    return (1.0 + std::erf(x / std::numbers::sqrt2)) / 2;
}

// Implementation of the Normal PDF function
double norm_pdf(double x)
{
    return std::exp(-0.5 * x * x) / std::sqrt(2 * std::numbers::pi);
}

// calculate the root of the function using Bisection and BSM
double bisection_method(BSM &bs, double market_price, bool debug)
{
    double a = 0.0001;
    double b = 3.0;
    double c = (a + b) / 2;
    double epsilon = 1e-06;
    double tol = bs(c) - market_price;
    int max_iter = 1000;
    int iter = 0;

    // **Check if the root is even bracketed**
    double fa = bs(a) - market_price;
    double fb = bs(b) - market_price;

    if (debug)
        std::cout << "bs(a) " << bs(a) << "bs(b) " << bs(b) << std::endl;

    if (fa * fb > 0)
    {
        if (debug)
            std::cout << "Warning: Root is not in range! Returning best estimate.\n";
        return std::min(std::max(epsilon, a), b);
    }

    if (debug)
        std::cout << "Starting Bisection Method...\n";

    while (std::abs(tol) > epsilon && iter < max_iter)
    {
        tol = bs(c) - market_price;

        if (tol < 0)
        {
            a = c;
        }
        else
        {
            b = c;
        }

        c = (a + b) / 2;
        iter++;

        // **Optional debug logging**
        if (debug && iter % 10 == 0)
        {
            std::cout << "Iteration " << iter << ": IV estimate = " << c << " | Tolerance = " << tol << std::endl;
        }
    }

    if (debug)
        std::cout << "Ending Bisection Method after " << iter << " iterations.\n";

    // **Ensure valid output within IV range**
    c = std::max(epsilon, std::min(c, b));

    return c;
}