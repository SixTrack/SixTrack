#pragma once
#include <vector>
#include "consts.h"
#include "signal.h"

std::vector<std::complex<double>> cheb_window(const size_t, const double);
std::vector<std::complex<double>> hann_harm_window(const size_t, const double);
std::vector<std::complex<double>> taylor_window(const size_t, const double);
std::vector<std::complex<double>> no_window(const size_t, const double);

class WindowFunc {
  private:
    Signal window;
  public:
    double parameter = 1.0;
    char type = 'h';

    WindowFunc(): window() {}
    ~WindowFunc() {}
    
    void compute(const size_t);
    void compute(const size_t, const double, const char);
    std::complex<double> operator()(size_t, size_t) const;
    std::complex<double> operator()(size_t, size_t);
};
