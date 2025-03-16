#include "IFD.h"
#include <cmath>
#include <iostream>
#include <vector>

IFD::IFD(double spot, double strike, double rate, double expiry, double sigma, OptionType option_type)
    : _spot(spot), _strike(strike), _rate(rate), _expiry(expiry), _sigma(sigma), _option_type(option_type)
{

    estimate_number_of_nodes();
}

void IFD::set_option_type(OptionType option_type)
{
    _option_type = option_type;
}

void IFD::estimate_number_of_nodes()
{
    double epsilon = 0.0001;

    // Set spatial grid size
    N = 1000;
    dx = 1.0 / N;

    // Compute dt using the implicit finite difference condition
    if (epsilon > dx * dx)
    {
        dt = epsilon - dx * dx; // Using equation (7.26)
    }
    else
    {
        std::cerr << "Error: dx^2 exceeds epsilon. Adjusting dt to a stable value." << std::endl;
        dt = dx * dx / 2.0; // Fallback stability condition
    }

    // Compute the number of time steps
    n = static_cast<int>(_expiry / dt);

    // Print computed values
    std::cout << "Implicit Finite Difference Grid:\n";
    std::cout << " dt: " << dt << "\n";
    std::cout << " n: " << n << "\n";
    std::cout << " dx: " << dx << "\n";
    std::cout << " N: " << N << "\n";
    std::cout << " dx^2 + dt = " << (dx * dx + dt) << std::endl;
}

double IFD::get_option_price()
{
    // compute  coefficients
    double nu = _rate - 0.5 * _sigma * _sigma;
    double A = -0.5 * dt * ((_sigma * _sigma) / (dx * dx) + nu / dx);
    double B = 1 + dt * (_sigma * _sigma) / (dx * dx);
    double C = -0.5 * dt * ((_sigma * _sigma) / (dx * dx) - nu / dx);

    double log_spot = std::log(_spot);

    // compute payoff at maturity
    std::vector<std::vector<double>> option_tree(n + 1, std::vector<double>(2 * N + 1, 0.0));
    for (int j = 0; j < 2 * N + 1; ++j)
    {
        double S = std::exp(log_spot + (j - N) * dx);
        if (_option_type == OptionType::EuropeanCall)
            option_tree[n][j] = std::max(S - _strike, 0.0);
        else
            option_tree[n][j] = std::max(_strike - S, 0.0);
    }

    // Boundary condition
    double lambda_U, lambda_L;
    if (_option_type == OptionType::EuropeanCall)
    {
        lambda_U = std::exp(N * dx + log_spot) * (std::exp(dx) - 1); // Upper boundary for Call
        lambda_L = 0.0;                                              // Lower boundary for Call
    }
    else // European Put
    {
        lambda_U = 0.0;                                               // Upper boundary for Put
        lambda_L = std::exp(-N * dx + log_spot) * (1 - std::exp(dx)); // Lower boundary for Put
    }

    // Backward time stepping
    for (int i = n - 1; i >= 0; i--)
    {
        solve_tridiagonal_matrix_gauss_elimination(option_tree[i], option_tree[i + 1], A, B, C, lambda_U, lambda_L);
    }

    // Return option price at S0 (center point at t=0)
    return option_tree[0][N];
}

void IFD::solve_tridiagonal_matrix(std::vector<double> &x,       // time i
                                   const std::vector<double> &y, // time i+1
                                   double A, double B, double C, double lambda_U, double lambda_L)
{

    int size = x.size();

    // Initialize vectors C and D
    std::vector<double> C_vector(size, 0.0);
    std::vector<double> D_vector(size, 0.0);

    // Apply upper boundary condition
    C_vector[0] = lambda_U;
    D_vector[0] = 1.0;

    // Forward elimination
    for (int i = 1; i < size - 1; i++)
    {
        double denom = B + A * D_vector[i - 1];
        C_vector[i] = (y[i] - A * C_vector[i - 1]) / denom;
        D_vector[i] = -C / denom;
    }

    // Apply lower boundary condition
    x[size - 1] = (lambda_L - A * C_vector[size - 2]) / (B + A * D_vector[size - 2]);

    // Backward substitution
    for (int i = size - 2; i >= 0; i--)
    {
        x[i] = C_vector[i] + D_vector[i] * x[i + 1];
    }
}

void IFD::solve_tridiagonal_matrix_gauss_elimination(std::vector<double> &x,       // time i
                                                     const std::vector<double> &y, // time i+1
                                                     double A, double B, double C, double lambda_U, double lambda_L)
{
    int size = x.size();

    // Modified coefficient vectors
    std::vector<double> alpha(size, 0.0); // Stores modified super-diagonal coefficients
    std::vector<double> beta(size, 0.0);  // Stores modified RHS values

    double epsilon = 1e-10; // To avoid division by zero errors

    // **Step 1: Forward elimination**
    alpha[0] = C / B;
    beta[0] = y[0] / B;

    for (int i = 1; i < size - 1; i++)
    {
        double denom = B - A * alpha[i - 1]; // Ensure denominator stability
        if (fabs(denom) < epsilon)
        {
            std::cerr << "Numerical instability detected at i = " << i << " (denom too small: " << denom << ")\n";
            return;
        }
        alpha[i] = C / denom;
        beta[i] = (y[i] - A * beta[i - 1]) / denom;
    }

    // **Step 2: Apply lower boundary condition**
    double denom = B - A * alpha[size - 2];
    if (fabs(denom) < epsilon)
    {
        std::cerr << "Numerical instability detected at boundary condition (denom too small: " << denom << ")\n";
        return;
    }
    x[size - 1] = (lambda_L - A * beta[size - 2]) / denom;

    // **Step 3: Backward substitution**
    for (int i = size - 2; i >= 0; i--)
    {
        x[i] = beta[i] - alpha[i] * x[i + 1];
    }
}

void IFD::estimate_grid_from_bs_price(double bs_price)
{

    // start first estimation with n = 1000
    n = 1000;
    N = n;
    dt = _expiry / n;
    dx = _sigma * std::sqrt(3 * dt);

    std::cout << "Estimation start with n: " << n << std::endl;
    std::cout << " dt: " << dt << "\n";
    std::cout << " dx: " << dx << "\n";
    std::cout << " N: " << N << "\n";

    double ifd_price = get_option_price();
    double epsilon = std::fabs(bs_price - ifd_price); // Ensure absolute difference
    int counter = 0;
    int max_iterations = 50; // Prevent infinite loops

    std::cout << "IFD price: " << ifd_price << std::endl;
    std::cout << "BSM price: " << bs_price << std::endl;

    while (epsilon > 0.001 && counter < max_iterations)
    {
        n *= 2; // Increase grid resolution
        N = n;
        dt = _expiry / n;
        dx = _sigma * std::sqrt(3 * dt);

        ifd_price = get_option_price();
        epsilon = std::fabs(bs_price - ifd_price); // Correct absolute error calculation
        counter++;

        std::cout << "Iteration " << counter << " - IFD price: " << ifd_price << ", Error: " << epsilon << std::endl;
    }

    std::cout << "Estimation end with n: " << n << std::endl;
    std::cout << "Final dt: " << dt << "\n";
    std::cout << "Final dx: " << dx << "\n";
    std::cout << "Final N: " << N << "\n";
    std::cout << "Total Iterations: " << counter << std::endl;

    if (counter >= max_iterations)
    {
        std::cerr << "Warning: Reached maximum iterations, convergence may not be achieved!" << std::endl;
    }
}
