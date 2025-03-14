#include "EFD.h"

/*
    Consider S0 = 100, K = 100, T = 1 year, σ = 20%, r = 6%, δ = 2%.
    Calculate and report the price for the European Call and Put using explicit,
    implicit, and Crank-Nicolson methods, while including the number of steps
    that you calculated in the previous point (part d).
*/

int main()
{

    double S0 = 100;
    double K = 100;
    double T = 1.0;
    double sigma = 0.2;
    double r = 0.06;
    double q = 0.02;

    EFD efd1(S0, K, r, T, sigma, OptionType::EuropeanCall);
}