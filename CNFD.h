#ifndef CNFD_H
#define CNFD_H

#include "OptionType.h"
#include <vector>

class CNFD
{
public:
    CNFD(double spot, double strike, double rate, double expiry, double sigma, OptionType option_type);
    double get_option_price();
    void set_option_type(OptionType option_type);
    void solve_tridiagonal_matrix(std::vector<double> &x, const std::vector<double> &y, double A, double B, double C, double lambda_U, double lambda_L);
    void solve_tridiagonal_matrix_gauss_elimination(std::vector<double> &x, const std::vector<double> &y, double A, double B, double C, double lambda_U, double lambda_L);
    void estimate_grid_from_bs_price(double bs_price);

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

#endif // CNFD_H
