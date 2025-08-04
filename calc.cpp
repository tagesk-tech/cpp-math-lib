#include <eigen3/Eigen/Dense>
#include <functional>
#include <iostream>
#include <cmath>        // for std::pow
#include <vector>
using namespace Eigen;

// f: Rⁿ → R (scalar)
// using Eigen::VectorXd; using Eigen::MatrixXd;

// compute gradient via central differences
VectorXd compute_gradient(const std::function<double(const VectorXd&)>& f,
                          const VectorXd& x, double h = 1e-6)
{
  int n = x.size();
  VectorXd grad(n);
  for(int i=0;i<n;++i){
    VectorXd xp = x, xm = x;
    xp[i] += h;  xm[i] -= h;
    grad[i] = (f(xp) - f(xm)) / (2*h);
  }
  return grad;
}

// compute Hessian via central differences
MatrixXd compute_hessian(const std::function<double(const VectorXd&)>& f,
                         const VectorXd& x, double h = 1e-6)
{
  int n = x.size();
  MatrixXd H(n,n);
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
      VectorXd xpp = x, xpm = x, xmp = x, xmm = x;
      xpp[i] += h; xpp[j] += h;
      xpm[i] += h; xpm[j] -= h;
      xmp[i] -= h; xmp[j] += h;
      xmm[i] -= h; xmm[j] -= h;
      H(i,j) = ( f(xpp) - f(xpm) - f(xmp) + f(xmm) ) / (4*h*h);
    }
  }
  return H;
}

VectorXd newtons_method(const std::function<double(const VectorXd&)>& f,
                        VectorXd x0, double tol = 1e-6, int max_iter = 100)
{
  VectorXd x = x0;
  for(int iter=0; iter<max_iter; ++iter){
    VectorXd g = compute_gradient(f, x);
    MatrixXd H = compute_hessian(f, x);
    // for SPD H:
    VectorXd dx = H.ldlt().solve(-g); //since the Hessian is symmetric you can quickly solve using LDL^T factorisation. 
    x += dx;
    if(dx.norm() < tol) break;
  }
  return x;
}


int main(){
    using Vec = Eigen::VectorXd;
    // f(x,y) = (x−3)² + (y+1)²  has ∇f=0 at (3,−1)
    auto f = [](const Vec& x){
      return std::pow(x[0] - 3.0, 2) //power function from the cmath library
             + std::pow(x[1] + 1.0, 2); 
           + std::pow(x[1] + 1.0, 2);
    };

    Vec x0(2);
    x0 << 0.0, 0.0;               // start at (0,0)
    Vec root = newtons_method(f, x0);

    std::cout << "Converged to: [" 
              << root[0] << ", " 
              << root[1] << "]\n";
    return 0;
}
