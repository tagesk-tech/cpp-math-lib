
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