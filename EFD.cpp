#include "EFD.h"
#include <cmath>
#include <iostream>
#include <vector>

EFD::EFD(double spot, double strike, double rate, double expiry, double sigma, OptionType option_type)
    : _spot(spot), _strike(strike), _rate(rate), _expiry(expiry), _sigma(sigma), _option_type(option_type)
{

    estimate_number_of_nodes();
}

void EFD::estimate_number_of_nodes()
{
    double epsilon = 0.0001;

    dt = epsilon / (1 + 3 * _sigma * _sigma);
    n = static_cast<int>(_expiry / dt); // Round to nearest integer
    N = n;
    dx = _sigma * std::sqrt(3 * dt);

    std::cout << "Calculated Numbers are: " << std::endl;
    std::cout << " dt: " << dt << std::endl;
    std::cout << " n: " << n << std::endl;
    std::cout << " N: " << N << std::endl;
    std::cout << " dx: " << dx << std::endl;
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
            option_tree[i][j] = A * option_tree[i + 1][j + 1] + B * option_tree[i + 1][j] + C * option_tree[i + 1][j - 1];
        }

        // apply Boundary Conditions
        option_tree[i][0] = option_tree[i][1];                                             // Lower boundary (approximation)
        option_tree[i][2 * N] = 2 * option_tree[i][2 * N - 1] - option_tree[i][2 * N - 2]; // Upper boundary
    }

    // Return option price at S_0 (central point)
    return option_tree[0][N];
}