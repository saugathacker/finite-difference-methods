#include "EFD.h"
#include "IFD.h"
#include "CNFD.h"
#include "BSM.h"
#include "FunctionsUtil.h"
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
    // // parameters
    // double S0 = 100;
    // double K = 100;
    // double T = 1.0;
    // double sigma = 0.2;
    // double r = 0.06;
    // double q = 0.02;

    // // creating FD objects with those paramters
    // EFD efd1(S0, K, r, T, sigma, OptionType::EuropeanCall);

    // IFD ifd1(S0, K, r, T, sigma, OptionType::EuropeanCall);

    // CNFD cnfd1(S0, K, r, T, sigma, OptionType::EuropeanCall);

    // // getting call prices for all three FD methods
    // double call_price = efd1.get_option_price();

    // cout << "The EFD call price is " << call_price << endl;

    // double ifd_call_price = ifd1.get_option_price();

    // cout << "The IFD call price is " << ifd_call_price << endl;

    // double cnfd_call_price = cnfd1.get_option_price();

    // cout << "The CNFD call price is " << cnfd_call_price << endl;

    // // changing the option type to put
    // efd1.set_option_type(OptionType::EuropeanPut);
    // ifd1.set_option_type(OptionType::EuropeanPut);
    // cnfd1.set_option_type(OptionType::EuropeanPut);

    // // getting all the put prices from FD methods
    // double put_price = efd1.get_option_price();

    // cout << "The EFD put price is " << put_price << endl;

    // double ifd_put_price = ifd1.get_option_price();

    // cout << "The IFD put price is " << ifd_put_price << endl;

    // double cnfd_put_price = cnfd1.get_option_price();

    // cout << "The CNFD put price is " << cnfd_put_price << endl;

    // calculating BSM price for estimating grid size for FD

    // BSM bsm1(K, S0, T, r, OptionType::EuropeanCall, 0.0);
    // double bsm_call_price = bsm1(sigma);

    // cout << "The BSM call price is " << bsm_call_price << endl;

    // efd1.estimate_grid_from_bs_price(bsm_call_price);

    // ifd1.estimate_grid_from_bs_price(bsm_call_price);

    // // test cases for solving tridiagobal matrix
    // int n = 3;
    // std::vector<double> x(n, 0.0);     // Unknowns
    // std::vector<double> y = {1, 2, 3}; // RHS
    // double A = -1.0;
    // double B = 2.0;
    // double C = -1.0;
    // double lambda_U = 1.0;
    // double lambda_L = 3.0;

    // ifd1.solve_tridiagonal_matrix(x, y, A, B, C, lambda_U, lambda_L);

    // std::cout << "\nSolution vector X:\n";
    // for (int i = 0; i < n; i++)
    // {
    //     std::cout << "x(" << i + 1 << ") = " << x[i] << std::endl;
    // }

    // ifd1.solve_tridiagonal_matrix_gauss_elimination(x, y, A, B, C, lambda_U, lambda_L);
    // std::cout << "\nSolution vector X:\n";
    // for (int i = 0; i < n; i++)
    // {
    //     std::cout << "x(" << i + 1 << ") = " << x[i] << std::endl;
    // }

    // int n2 = 5;
    // std::vector<double> x2(n2, 0.0);            // Unknowns
    // std::vector<double> y2 = {5, 5, 10, 10, 5}; // RHS
    // double A2 = -1.0;                           // Sub-diagonal
    // double B2 = 4.0;                            // Main diagonal
    // double C2 = -1.0;                           // Super-diagonal
    // double lambda_U2 = 5.0;                     // Upper boundary condition
    // double lambda_L2 = 5.0;                     // Lower boundary condition

    // ifd1.solve_tridiagonal_matrix_gauss_elimination(x2, y2, A2, B2, C2, lambda_U2, lambda_L2);

    // std::cout << "\nSolution vector X:\n";
    // for (int i = 0; i < n2; i++)
    // {
    //     std::cout << "x(" << i + 1 << ") = " << x2[i] << std::endl;
    // }

    // double delta, gamma, theta, vega;
    // efd1.compute_greeks(delta, gamma, theta, vega);

    // std::cout << "Delta: " << delta << "\n";
    // std::cout << "Gamma: " << gamma << "\n";
    // std::cout << "Theta: " << theta << "\n";
    // std::cout << "Vega: " << vega << "\n";

    // from the downloaded data calculating the IV and using that to get FD prices
    // FunctionsUtil fu;
    // std::unordered_map<std::string, std::unique_ptr<Ticker>> tickers_data = fu.read_csv_into_ticker_object("NVDA_options_data.csv");
    // fu.perform_operations(tickers_data);
    return 0;
}