#include <eigen3/Eigen/Dense>
#include <functional>
#include <iostream>
#include <algorithm>

using Vec = Eigen::VectorXd;
using Mat = Eigen::MatrixXd;

// Finite-difference gradient of scalar f
Vec compute_gradient(const std::function<double(const Vec&)>& f,
                     const Vec& x, double h = 1e-6) {
    int n = x.size(); Vec g(n);
    for(int i = 0; i < n; ++i) {
        Vec xp = x, xm = x;
        xp(i) += h; xm(i) -= h;
        g(i) = (f(xp) - f(xm)) / (2*h);
    }
    return g;
}

// Finite-difference Hessian of scalar f
Mat compute_hessian(const std::function<double(const Vec&)>& f,
                    const Vec& x, double h = 1e-4) {
    int n = x.size(); Mat H(n,n);
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            Vec xpp = x, xpm = x, xmp = x, xmm = x;
            xpp(i)+=h; xpp(j)+=h;
            xpm(i)+=h; xpm(j)-=h;
            xmp(i)-=h; xmp(j)+=h;
            xmm(i)-=h; xmm(j)-=h;
            H(i,j) = ( f(xpp) - f(xpm) - f(xmp) + f(xmm) ) / (4*h*h);
        }
    }
    return H;
}

// Finite-difference Jacobian of vector-valued g(x)
Mat compute_jacobian(const std::function<Vec(const Vec&)>& g,
                     const Vec& x, double h = 1e-6) {
    Vec y = g(x);
    int m = y.size(), n = x.size();
    Mat J(m,n);
    for(int i = 0; i < n; ++i) {
        Vec xp = x, xm = x;
        xp(i) += h; xm(i) -= h;
        J.col(i) = (g(xp) - g(xm)) / (2*h);
    }
    return J;
}

// Problem definition
struct NLP {
    std::function<double(const Vec&)> f; // objective
    std::function<Vec(const Vec&)> g;    // inequality constraints g(x) <= 0
    std::function<Vec(const Vec&)> h;    // equality constraints h(x) = 0
};

// Primal-dual Newton (KKT) solver
bool solve_kkt(const NLP& P,
               Vec& x, Vec& lambda, Vec& nu,
               double tol = 1e-6, int max_iter = 50) {
    int n = x.size(), m = lambda.size(), p = nu.size();
    double rho = 1.0;
    for(int iter = 0; iter < max_iter; ++iter) {
        // Residuals
        Vec df = compute_gradient(P.f, x);
        Mat G = compute_jacobian(P.g, x);
        Mat Hc = compute_jacobian(P.h, x);
        Vec gi = P.g(x);
        Vec hi = P.h(x);

        // Hessian of Lagrangian
        Mat H = compute_hessian(P.f, x);
        for(int i = 0; i < m; ++i)
            H += lambda(i) * compute_hessian(
                  [&](const Vec& v){ return P.g(v)(i); }, x);
        for(int j = 0; j < p; ++j)
            H += nu(j) * compute_hessian(
                  [&](const Vec& v){ return P.h(v)(j); }, x);

        // Build KKT matrix and residual F
        int sz = n + m + p;
        Mat K = Mat::Zero(sz, sz);
        Vec F = Vec::Zero(sz);

        // Top blocks
        K.block(0,0,n,n) = H;
        K.block(0,n,n,m)     = G.transpose();
        K.block(0,n+m,n,p)   = Hc.transpose();
        F.segment(0,n)       = df + G.transpose()*lambda + Hc.transpose()*nu;

        // Complementarity row
        K.block(n,0,m,n)     = lambda.asDiagonal() * G;
        K.block(n,n,m,m)     = gi.asDiagonal();
        F.segment(n,m)       = lambda.array() * gi.array();

        // Equality row
        K.block(n+m,0,p,n)   = Hc;
        F.segment(n+m,p)     = hi;

        // Solve for search direction
        Vec dz = K.fullPivLu().solve(-F);
        Vec dx     = dz.segment(0,n);
        Vec dlambda= dz.segment(n,m);
        Vec dnu    = dz.segment(n+m,p);

        // Fraction-to-boundary for lambda>=0
        double alpha = 1.0;
        for(int i = 0; i < m; ++i)
            if(dlambda(i) < 0)
                alpha = std::min(alpha, -lambda(i)/dlambda(i)*0.995);

        // Merit function
        auto phi = [&](const Vec& xx){
            Vec gg = P.g(xx);
            Vec hh = P.h(xx);
            return P.f(xx)
                 + 0.5*rho*(hh.squaredNorm() + gg.cwiseMax(0.0).squaredNorm());
        };
        double phi0 = phi(x);
        Vec grad_phi = compute_gradient(phi, x);

        // Backtracking line search
        double c = 1e-4;
        while(phi(x + alpha*dx)
              > phi0 + c*alpha*grad_phi.dot(dx))
            alpha *= 0.5;

        // Update primal and dual
        x      += alpha * dx;
        lambda += alpha * dlambda;
        nu     += alpha * dnu;

        // Check convergence
        if(F.norm() < tol)
            return true;
    }
    return false;
}

int main() {
    // Example: minimize (x-3)^2 s.t. x<=2
    NLP P;
    P.f = [](const Vec& v){ return std::pow(v[0] - 3.0, 2); };
    P.g = [](const Vec& v){ return Vec::Constant(1, v[0] - 2.0); };
    P.h = [](const Vec&){ return Vec(); };

    Vec x(1); x << 0.0;
    Vec lambda = Vec::Ones(1);
    Vec nu     = Vec();

    bool ok = solve_kkt(P, x, lambda, nu);
    if(ok)
        std::cout << "Solution x=" << x.transpose()
                  << ", lambda=" << lambda.transpose() << std::endl;
    else
        std::cout << "Did not converge" << std::endl;

    return 0;
}
