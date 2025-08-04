// This file is part of the C++ library for basic calculus operations. 
// like gradients, hessians, and other calculus related operations.  

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

// example on how the function and gradient and hessian should look like


std::vector<double> compute_gradient(
    const std::function<double(const std::vector<double>&)>& f, // the function to compute the gradient of
    int n, // the number of variables in the function
    const std::vector<double>& x0, // the point at which to compute the gradient
    double h = 1e-6 // the step size for finite difference approximation
) {
    std::vector<double> gradient(n);
    for (int i = 0; i < n; ++i) {
        std::vector<double> x_plus = x0;
        std::vector<double> x_minus = x0;
        x_plus[i] += h;
        x_minus[i] -= h;
        gradient[i] = (f(x_plus) - f(x_minus)) / (2 * h);
    }
    return gradient;   
}

// calculates teh hessian matrix of a function at a given point using finite difference approximation
std::vector<std::vector<double>> compute_hessian(
    const std::function<double(const std::vector<double>&)>& f, // the function to compute the hessian of
    int n, // the number of variables in the function
    const std::vector<double>& x0, // the point at which to compute the hessian
    double h = 1e-6 // the step size for finite difference approximation
) {
    std::vector<std::vector<double>> hessian(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::vector<double> x_plus_i = x0;
            std::vector<double> x_minus_i = x0;
            std::vector<double> x_plus_j = x0;
            std::vector<double> x_minus_j = x0;
            x_plus_i[i] += h;
            x_minus_i[i] -= h;
            x_plus_j[j] += h;
            x_minus_j[j] -= h;

            double f_plus_plus = f(x_plus_i);
            double f_plus_minus = f(x_plus_j);
            double f_minus_plus = f(x_minus_i);
            double f_minus_minus = f(x_minus_j);

            hessian[i][j] = (f_plus_plus - f_plus_minus - f_minus_plus + f_minus_minus) / (4 * h * h);
        }
    }
    return hessian;
}
