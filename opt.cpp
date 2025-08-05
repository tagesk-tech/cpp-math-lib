
#pragma once
#include "calc.h" // in other .cpp files, you can include this header to use the functions defined in calc.h
#include <iostream>
#include <vector>
#include <stdexcept>    // for std::runtime_error
#include <algorithm>    // for std::swap
#include <eigen3/Eigen/Dense> // remember eigen3
using namespace Eigen;
using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;


// Our goal in this case is to solve non liniar optimization problems using the Newton's method


struct NLProblem {
  std::function<double(const Vec&)>      f;  // objective
  std::function<Vec  (const Vec&)>      g;  // g(x): Rⁿ→Rᵐ (inequalities)
  std::function<Vec  (const Vec&)>      h;  // h(x): Rⁿ→Rᵖ (equalities)

  NLProblem( decltype(f) _f,
             decltype(g) _g,
             decltype(h) _h )
    : f(std::move(_f))
    , g(std::move(_g))
    , h(std::move(_h))
  {}

    // we will start out with a function that is convex and assume it is correct in this case

    Vec solve(const Vec& x0, double tol = 1e-6, int max_iter = 100) const {

};

};




int main(){
    return 0; // just a placeholder for now
}

