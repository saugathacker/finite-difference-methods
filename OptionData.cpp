#include "OptionData.h"
#include <iostream>
#include <cmath>
#include "util.h"
#include "BSM.h"
#include "EFD.h"
#include "IFD.h"
#include "CNFD.h"

// implementation of calculate iv and greeks
void OptionData::calculate_iv(double spotPrice, double interestRate)
{
    // Skip calculation if lastPrice, bid, and ask are all zero
    if ((lastPrice <= 0.0 || std::isnan(lastPrice)) &&
        (bid <= 0.0 || std::isnan(bid)) &&
        (ask <= 0.0 || std::isnan(ask)))
    {
        return; // Do not calculate IV if no valid market data exists
    }

    // Use mid-price if bid/ask exist, otherwise use last price
    double market_price = (bid > 0 && ask > 0) ? (bid + ask) / 2 : lastPrice;

    // Create Black-Scholes model instance
    BSM bs(strike, spotPrice, timeToMaturity, interestRate,
           (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));

    // std::cout << bs << std::endl;
    // std::cout << "Finding root with market price: " << market_price << std::endl;
    // Measure time for Bisection Method
    auto start_bisect = std::chrono::high_resolution_clock::now();
    bisectionImpliedVol = bisection_method(bs, market_price, false);
    auto end_bisect = std::chrono::high_resolution_clock::now();
    bisectionTime = std::chrono::duration<double, std::milli>(end_bisect - start_bisect).count();
}

void OptionData::calculate_fd_price(double spot, double rate)
{
    EFD efd(spot, strike, rate, timeToMaturity, bisectionImpliedVol,
            (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));
    efd_price = efd.get_option_price();

    IFD ifd(spot, strike, rate, timeToMaturity, bisectionImpliedVol,
            (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));
    ifd_price = ifd.get_option_price();

    CNFD cnfd(spot, strike, rate, timeToMaturity, bisectionImpliedVol,
              (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));
    cnfd_price = cnfd.get_option_price();
}

void OptionData::calculate_greeks_efd(double spot, double rate)
{
    EFD efd(spot, strike, rate, timeToMaturity, bisectionImpliedVol,
            (optionType == "Call" ? OptionType::EuropeanCall : OptionType::EuropeanPut));
    efd_price = efd.get_option_price();

    efd.compute_greeks(efd_delta, efd_gamma, efd_theta, efd_vega);
}