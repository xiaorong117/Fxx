#include <badiff.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

double function(double x, double y) {
  return x * x * std::exp(y) + std::sin(x * y);
}

double relative_or_absolute(double actual, double expected) {
  return std::abs(actual - expected) / std::max(std::abs(expected), 1.0e-300);
}

}  // namespace

int main() {
  const double xv = 0.7;
  const double yv = -0.2;
  fadbad::B<double> x = xv, y = yv;
  fadbad::B<double> f = x * x * exp(y) + sin(x * y);
  f.diff(0, 1);
  const double ad_x = x.d(0);
  const double ad_y = y.d(0);
  const double exact_x = 2.0 * xv * std::exp(yv) + yv * std::cos(xv * yv);
  const double exact_y = xv * xv * std::exp(yv) + xv * std::cos(xv * yv);
  const double hx = std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, std::abs(xv));
  const double hy = std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, std::abs(yv));
  const double fd_x = (function(xv + hx, yv) - function(xv - hx, yv)) / (2.0 * hx);
  const double fd_y = (function(xv, yv + hy) - function(xv, yv - hy)) / (2.0 * hy);
  std::cout << std::scientific << std::setprecision(17)
            << "dfdx ad=" << ad_x << " analytic=" << exact_x << " fd=" << fd_x << '\n'
            << "dfdy ad=" << ad_y << " analytic=" << exact_y << " fd=" << fd_y << '\n';
  if (relative_or_absolute(ad_x, exact_x) > 1e-9 ||
      relative_or_absolute(ad_y, exact_y) > 1e-9 ||
      relative_or_absolute(ad_x, fd_x) > 1e-9 ||
      relative_or_absolute(ad_y, fd_y) > 1e-9) {
    throw std::runtime_error("FADBAD++ smoke derivative tolerance exceeded");
  }
  return 0;
}
