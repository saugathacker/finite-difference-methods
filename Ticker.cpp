#include "Ticker.h"
#include <fstream>
#include <iostream>

// Constructor
Ticker::Ticker(const std::string &name, double spot, double rate)
    : tickerName(name), spotPrice(spot), interestRate(rate) {}

// add new option data to the existing Ticker object
void Ticker::addOptionData(std::unique_ptr<OptionData> option)
{
    options.push_back(std::move(option));
}

// find an option based on strike, expiration, and type
OptionData *Ticker::findOption(double strike, const std::string &expiration, const std::string &optionType) const
{
    for (const auto &option : options)
    {
        if (option->strike == strike && option->expiration == expiration && option->optionType == optionType)
        {
            return option.get(); // return raw pointer (safe since unique_ptr manages memory)
        }
    }
    return nullptr; // return nullptr if no matching option is found
}

void Ticker::calculate_implied_vol_from_bs()
{
    for (auto &option : options)
    {
        if (option->strike >= 100 && option->strike <= 140)
        {
            option->calculate_iv(spotPrice, interestRate); // each option calculates its IV
        }
    }
}

void Ticker::calculate_fd_prices()
{
    for (auto &option : options)
    {
        if (option->strike >= 100 && option->strike <= 140)
        {
            option->calculate_fd_price(spotPrice, interestRate);
        }
    }
}

void Ticker::calculate_greeks_efd()
{
    for (auto &option : options)
    {
        if (option->strike >= 100 && option->strike <= 140)
        {
            option->calculate_greeks_efd(spotPrice, interestRate);
        }
    }
}

// implentaion of write to csv all the option data (observed and calculated) of this Ticker
void Ticker::write_to_csv(const std::string &filename) const
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // **Write CSV Header**
    file << "Ticker,Expiration,TimeToMaturity,Strike,OptionType,LastPrice,"
         << "Bid,Ask,ImpliedVolatility,BisectionIV,BisectionTime,EFDprice,IFDprice,CNFDprice,EFDdelta,EFDgamma,EFDtheta,EFDvega,InTheMoney\n";

    // **Write Option Data**
    for (const auto &option : options)
    {
        file << tickerName << "," // Add ticker symbol
             << option->expiration << ","
             << option->timeToMaturity << ","
             << option->strike << ","
             << option->optionType << ","
             << option->lastPrice << ","
             << option->bid << ","
             << option->ask << ","
             << option->impliedVolatility << ","
             << option->bisectionImpliedVol << ","
             << option->bisectionTime << ","
             << option->efd_price << ","
             << option->ifd_price << ","
             << option->cnfd_price << ","
             << option->efd_delta << ","
             << option->efd_gamma << ","
             << option->efd_theta << ","
             << option->efd_vega << ","
             << (option->inTheMoney ? "True" : "False") << "\n";
    }

    file.close();
    std::cout << "CSV file written successfully: " << filename << std::endl;
}