#include "EFD.h"
#include "BSM.h"
#include <iostream>

using namespace std;

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

    double call_price = efd1.get_option_price();

    cout << "The EFD call price is " << call_price << endl;

    BSM bsm1(K, S0, T, r, OptionType::EuropeanCall, 0.0);
    double bsm_call_price = bsm1(sigma);

    cout << "The BSM call price is " << bsm_call_price << endl;

    return 0;
}