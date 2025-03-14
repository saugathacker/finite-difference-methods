#pragma once

#include "OptionType.h"

class EFD
{
public:
    EFD(double spot, double strike, double rate, double expiry, double sigma, OptionType option_type);
    double get_option_price();

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

    void estimate_number_of_nodes();
};
