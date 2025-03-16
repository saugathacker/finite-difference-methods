#pragma once

#include "OptionType.h"
#include <vector>

class EFD
{
public:
    EFD(double spot, double strike, double rate, double expiry, double sigma, OptionType option_type);
    double get_option_price();
    void estimate_grid_from_bs_price(double bs_price);
    void compute_greeks(double &delta, double &gamma, double &theta, double &vega);
    void set_option_type(OptionType option_type);

private:
    double _spot;
    double _strike;
    double _rate;
    double _expiry;
    double _sigma;

    int n;
    int N;
    double dt;
    double dx;

    OptionType _option_type;
    std::vector<std::vector<double>> _option_tree;

    void estimate_number_of_nodes();
};
