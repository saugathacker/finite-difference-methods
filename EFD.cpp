#include "EFD.h"
#include <cmath>
#include <iostream>
#include <vector>

EFD::EFD(double spot, double strike, double rate, double expiry, double sigma, OptionType option_type)
    : _spot(spot), _strike(strike), _rate(rate), _expiry(expiry), _sigma(sigma), _option_type(option_type)
{

    estimate_number_of_nodes();
}

void EFD::set_option_type(OptionType option_type)
{
    _option_type = option_type;
}

void EFD::estimate_number_of_nodes()
{
    double epsilon = 0.0001;

    dt = epsilon / (1 + 3 * _sigma * _sigma);
    n = static_cast<int>(_expiry / dt + 0.5); // Round to nearest integer
    N = n;
    dx = _sigma * std::sqrt(3 * dt);

    std::cout << "Explicit Finite Difference Grid:\n";
    std::cout << " dt: " << dt << "\n";
    std::cout << " n: " << n << "\n";
    std::cout << " dx: " << dx << "\n";
    std::cout << " N: " << N << "\n";
    std::cout << " sigma^2 * 3 * dt + dt = " << (3.0 * _sigma * _sigma * dt + dt) << std::endl;
}

double EFD::get_option_price()
{
    // compute  coefficients
    double nu = _rate - 0.5 * _sigma * _sigma;
    double A = 0.5 * dt * ((_sigma * _sigma) / (dx * dx) + nu / dx);
    double B = 1 - dt * (_sigma * _sigma) / (dx * dx);
    double C = 0.5 * dt * ((_sigma * _sigma) / (dx * dx) - nu / dx);

    double log_spot = std::log(_spot);

    // compute payoff at maturity
    std::vector<std::vector<double>> option_tree(n, std::vector<double>(2 * N + 1, 0.0));
    for (int j = 0; j < 2 * N + 1; ++j)
    {
        double S = std::exp(log_spot + (j - N) * dx);
        if (_option_type == OptionType::EuropeanCall)
            option_tree[n - 1][j] = std::max(S - _strike, 0.0);
        else
            option_tree[n - 1][j] = std::max(_strike - S, 0.0);
    }

    // backward induction
    for (int i = n - 2; i >= 0; --i)
    {
        for (int j = 1; j < 2 * N; ++j)
        { // Avoid boundaries for now
            option_tree[i][j] = (1 / (1 + (_rate * dt))) * (A * option_tree[i + 1][j + 1] +
                                                            B * option_tree[i + 1][j] + C * option_tree[i + 1][j - 1]);
        }

        // apply Boundary Conditions
        option_tree[i][0] = option_tree[i][1];                                             // Lower boundary (approximation)
        option_tree[i][2 * N] = 2 * option_tree[i][2 * N - 1] - option_tree[i][2 * N - 2]; // Upper boundary
    }

    _option_tree = option_tree;

    // Return option price at S_0 (central point)
    return option_tree[0][N];
}

void EFD::estimate_grid_from_bs_price(double bs_price)
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

    double efd_price = get_option_price();
    double epsilon = std::fabs(bs_price - efd_price); // Ensure absolute difference
    int counter = 0;
    int max_iterations = 50; // Prevent infinite loops

    std::cout << "EFD price: " << efd_price << std::endl;
    std::cout << "BSM price: " << bs_price << std::endl;

    while (epsilon > 0.001 && counter < max_iterations)
    {
        n *= 2; // Increase grid resolution
        N = n;
        dt = _expiry / n;
        dx = _sigma * std::sqrt(3 * dt);

        efd_price = get_option_price();
        epsilon = std::fabs(bs_price - efd_price); // Correct absolute error calculation
        counter++;

        std::cout << "Iteration " << counter << " - EFD price: " << efd_price << ", Error: " << epsilon << std::endl;
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

void EFD::compute_greeks(double &delta, double &gamma, double &theta, double &vega)
{
    double original_price = get_option_price();
    int t0 = 0; // Current time step

    // Compute finite step size in actual price terms
    double S_up = std::exp(dx + std::log(_spot));
    double S_mid = _spot;
    double S_down = std::exp(-dx + std::log(_spot));

    // **Delta Calculation**
    delta = (_option_tree[t0][N + 1] - _option_tree[t0][N - 1]) / (S_up - S_down);

    // **Gamma Calculation**
    gamma = (_option_tree[t0][N + 1] - 2 * _option_tree[t0][N] + _option_tree[t0][N - 1]) / ((S_up - S_mid) * (S_mid - S_down));

    // **Theta Calculation** (Forward Difference)
    theta = (_option_tree[t0 + 1][N] - _option_tree[t0][N]) / dt;

    // **Vega Calculation** (Recomputing Option Price with Perturbed Sigma)
    double sigma_original = _sigma;
    _sigma += 0.0001;                       // Small increment in volatility
    double price_vega = get_option_price(); // Recalculate option price
    _sigma = sigma_original;                // Restore original sigma
    vega = (price_vega - original_price) / 0.0001;
}
