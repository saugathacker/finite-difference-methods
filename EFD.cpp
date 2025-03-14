#include "EFD.h"
#include <cmath>
#include <iostream>

EFD::EFD(double spot, double strike, double rate, double expiry, double sigma, OptionType option_type)
    : _spot(spot), _strike(strike), _rate(rate), _expiry(expiry), _sigma(sigma), _option_type(option_type)
{

    estimate_number_of_nodes();
}

void EFD::estimate_number_of_nodes()
{
    double epsilon = 0.0001;

    dt = epsilon / (1 + 3 * _sigma * _sigma);
    n = _expiry / dt;
    N = n;
    dx = _sigma * std::sqrt(3 * dt);

    std::cout << "Calculated Numbers are: " << std::endl;
    std::cout << " dt: " << dt << std::endl;
    std::cout << " n: " << n << std::endl;
    std::cout << " N: " << N << std::endl;
    std::cout << " dx: " << dx << std::endl;
}