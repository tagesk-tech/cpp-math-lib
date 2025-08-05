#include <iostream>
#include <vector>
#include <cstdio>        // for popen, fprintf
#include <eigen3/Eigen/Dense>
#include <functional>
#include <cmath>        // for std::pow
#include <stdexcept>    // for std::runtime_error
#include <algorithm>    // for std::swap



using namespace Eigen;
using Vec = VectorXd;
using Mat = MatrixXd;

// Finite-difference solver for a simple boundary-value problem
enum { PLOT_USE_PERSIST = 1 };
Vec finite_difference(int n, double y0, double yn, double start, double end) { // finite difference only works for non-liniar problems
    Mat A = Mat::Zero(n, n);
    Vec b = Vec::Zero(n);
    double h = (end - start) / (n + 1);

    // Build tridiagonal matrix A
    for (int i = 0; i < n; ++i) {
        A(i, i) = -2 - 4 * h * h;
        if (i > 0)     A(i, i - 1) = 1;
        if (i < n - 1) A(i, i + 1) = 1;
    }
    // Right-hand side with boundary conditions y(0)=y0, y(end)=yn
    b(0)     = -y0;
    b(n - 1) = -yn;

    // Solve A y = b
    return A.ldlt().solve(b);
}

Vec newton_system(
  const std::function<Vec(const Vec&)>& F,
  const std::function<Mat(const Vec&)>& J,
  Vec                                  x0,
  double                               tol      = 1e-6,
  int                                  max_iter = 100
){
  Vec x = x0;
  for(int iter=0; iter<max_iter; ++iter){
    Vec Fx = F(x);
    Mat Jx = J(x);
    Vec dx = Jx.ldlt().solve(-Fx);
    x += dx;
    if(dx.norm() < tol) break;
  }
  return x;
}

Vec F(const Vec& x, int start, int end, const Vec& bv, int n){
  Vec y = Vec::Zero(n);
  double h  = double(end-start)/(n+1), hh = h*h;
  // first interior node
  y(0) = bv(0)
       - (2+hh)*x(0)
       + x(1)
       +     hh * x(0)*x(0);
  // interior
  for(int i=1; i<n-1; ++i){
    y(i) = x(i-1)
         - (2+hh)*x(i)
         + x(i+1)
         +     hh * x(i)*x(i);
  }
  // last interior node
  y(n-1) = x(n-2)
         - (2+hh)*x(n-1)
         +   bv(1)
         +     hh * x(n-1)*x(n-1);
  return y;
}

Mat J(const Vec& x, int start, int end, const Vec& /*bv*/, int n){
  Mat J = Mat::Zero(n,n);
  double h  = double(end-start)/(n+1), hh = h*h;
  // first row
  J(0,0)   = -(2+hh) + 2*hh*x(0);
  if(n>1) J(0,1)   = 1;
  // interior
  for(int i=1; i<n-1; ++i){
    J(i,i-1) = 1;
    J(i,i)   = -(2+hh) + 2*hh*x(i);
    J(i,i+1) = 1;
  }
  // last row
  if(n>1){
    J(n-1,n-2) = 1;
    J(n-1,n-1) = -(2+hh) + 2*hh*x(n-1);
  }
  return J;
}

Vec finite_difference_NLP(int start, int end, Vec bv, int n){
  auto Ff = [&](const Vec& x){ return F(x,start,end,bv,n); };
  auto Jf = [&](const Vec& x){ return J(x,start,end,bv,n); };
  Vec x0 = Vec::Zero(n);
  x0(0)=bv(0); x0(n-1)=bv(1);
  return newton_system(Ff,Jf,x0);
}

Vec heat(){
    // This function is a placeholder for the heat equation solver.
    // You can implement a specific heat equation solver here if needed.
    return Vec();
}

// Plot (x,y) using gnuplot via a pipe
void plot_solution(const Vec& y, double start, double end) {
    int n = y.size();
    double h = (end - start) / (n + 1);
    FILE* gp = popen("gnuplot -persist", "w");
    if (!gp) {
        std::cerr << "Error: could not open pipe to gnuplot.\n";
        return;
    }

    fprintf(gp, "set title 'Finite Difference Solution'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "plot '-' with linespoints title 'y(x)'\n");

    for (int i = 0; i < n; ++i) {
        double x = start + (i + 1) * h;
        fprintf(gp, "%g %g\n", x, y[i]);
    }
    fprintf(gp, "e\n");
    fflush(gp);
    // leave gp open if you want the plot window to persist
}



int main() {
    int n      = 10;      // 10 interior points
    double start = 0.0;
    double end   = 1.0;

    // boundary values y(0)=1, y(1)=10
    Vec bv(2);
    bv << 1.0, 4.0;

    // solve the non-linear finite-difference system
    Vec solution = finite_difference_NLP(start, end, bv, n);

    // plot with gnuplot
    plot_solution(solution, start, end);

    return 0;
}
