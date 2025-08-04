#include <iostream>
#include <vector>
#include <initializer_list>
#include <sstream>
#include <string>
// polynomial.cpp


// målet är enkelt att skriva ett program som kan jobba med polynom och som framörallt kan hitta GCD och LCM av polynom
// splita det i chinese remainder theorem delar också.


// modulo class that has the modulo property to it in this case. :)
// målet är att ha ett kostnadseffektivt bibliotek för att jobba med modulo
#include <iostream>

template <int MOD>
class Mod {
    int v;
public:
    // construct from any integer, reduce into [0..MOD)
    explicit Mod(int x = 0) // lets it still be an int in this case and does not convert it into a Mod object in this case
      : v((x % MOD + MOD) % MOD) // here is just how c++ lets you say how each member gets it initial value
    {}

    // accessors
    int value() const { return v; }

    // operator+= and *= do the work then reduce
    Mod& operator+=(Mod const& o) { //Mod& means i am returning a reference to the same object you called this on
        v += o.v;
        if (v >= MOD) v -= MOD;
        return *this; //hand you back the original thing i called this on
    }
    
    Mod& operator*=(Mod const& o) {
        v = int((1LL * v * o.v) % MOD); //avoids overflow since i am stating basically compute this as a long.
        return *this;
    }

    // binary + and * in terms of the compound ones
    friend Mod operator+(Mod a, Mod const& b) {  // the friend operator lets you peek into the internal things of the class
        return a += b; 
    }
    friend Mod operator*(Mod a, Mod const& b) { 
        return a *= b; 
    }

    // (optional) stream-output
    friend std::ostream& operator<<(std::ostream& os, Mod const& m) {
        return os << m.v;
    }

    // 1) Make sure your Mod<MOD> supports == and !=
    friend bool operator==(Mod const& a, Mod const& b) { return a.v == b.v; }
    friend bool operator!=(Mod const& a, Mod const& b) { return a.v != b.v; }
};

template<typename Field>
class Polynomial {
    std::vector<Field> coeffs;  // coeffs[i] is the x^i coefficient // you seem to be able to work with the FFT here (absolut värt det för att skynda på processen)

public:
  // my constructor function in this case
  Polynomial(std::initializer_list<Field> init)
    : coeffs(init)
  {
    trim();  // remove any trailing zeros
  }
    
    // constructors, e.g. from an initializer_list<Field>
    // + operator: pad the shorter coeff vector, add termwise via Field::operator+=
    // * operator: convolution using Field::operator*= and operator+=
    // evaluation, differentiation, etc.

    void trim() { // remove trailing zeros
        while (!coeffs.empty() && coeffs.back() == Field(0)) {
            coeffs.pop_back();
        }
    }


    Polynomial& operator+=(Polynomial const& other) {
        // make sure the coeffs vector is large enough
        coeffs.resize(std::max(coeffs.size(), other.coeffs.size()));
        for (size_t i = 0; i < other.coeffs.size(); ++i) {
            coeffs[i] += other.coeffs[i]; // add the coefficients together
        }
        return *this; // return the current object
    }

    Polynomial& operator*=(Polynomial const& other) {
        std::vector<Field> result(coeffs.size() + other.coeffs.size() - 1); // result size is the sum of the sizes of a and b minus 1
        for (size_t i = 0; i < coeffs.size(); ++i) {
            for (size_t j = 0; j < other.coeffs.size(); ++j) {
                result[i + j] += coeffs[i] * other.coeffs[j]; // multiply the coefficients and add them to the result pretty standard
            }
        }
        coeffs = std::move(result); // move the result into coeffs
        return *this; // return the current object
    }

    Field eval(Field x) const { //evaluates the polynomial at a point x using the field
        Field y{0};
        for (int i = coeffs.size() - 1; i >= 0; --i) { // evaluate the polynomial at x
            y = y * x + coeffs[i]; // Horner's method
        }
        return y;    // returns the Field value in this case
    }
    template<typename F>
    static std::vector<F> all_field_elements() {
    std::vector<F> elems;
    F cur{0};
    do {
        elems.push_back(cur);
        cur = cur + F{1};
    } while (cur != F{0});
    return elems;
    }


    // 4) A brute‐force reduce() that peels off any linear factor x–α you find:
    // It could be interesting to implement a reduce() method that finds and removes linear factors of the polynomial. and then factors the polynomial accordingly.
    Polynomial& reduce() {
        auto elems = all_field_elements<Field>();
        for (auto const& a : elems) {
            if (eval(a) == Field{0}) {
                // we found a root a, so x–a is a factor
                // form the “x – a” polynomial:
                Polynomial<Field> lin{ Field{-a.value()}, Field{1} };
                // now do a naive polynomial division of *this by lin:
                std::vector<Field> q(coeffs.size()), r = coeffs;
                for (int i = int(r.size()) - 1; i >= 1; --i) {
                    q[i-1] = r[i];
                    for (int j = 0; j < 2; ++j)
                        r[i-1-j] += q[i-1] * lin.coeffs[1-j];
                }
                q.resize(r.size()-1);
                coeffs = std::move(q);
                trim();
                break;   // if you want just one step; loop to fully factor
            }
        }
        return *this;
    }


    
    friend std::vector<Field> operator+(Polynomial const& a, Polynomial const& b) { //we dont want to change the value for a or b in this case
        std::vector<Field> result(std::max(a.coeffs.size(), b.coeffs.size())); // create a result vector that is the size of the largest coeffs vector
        for (size_t i = 0; i < result.size(); ++i) {
            if (i < a.coeffs.size()) result[i] += a.coeffs[i]; // if i is less than the size of a.coeffs then add the value of a.coeffs[i] to result[i]
            if (i < b.coeffs.size()) result[i] += b.coeffs[i];
        }
        return result;
    }

    friend std::vector<Field> operator*(Polynomial const& a, Polynomial const& b){

        std::vector<Field> result(a.coeffs.size() + b.coeffs.size() - 1); // result size is the sum of the sizes of a and b minus 1
        for (size_t i = 0; i < a.coeffs.size(); ++i) {
            for (size_t j = 0; j < b.coeffs.size(); ++j) {
                result[i + j] += a.coeffs[i] * b.coeffs[j]; // multiply the coefficients and add them to the result pretty standard
            }
        }
        return result;
    }



};


int main() {

    using F5 = Mod<7>; // change your field here

    // 1) A reducible polynomial: x² + 1 over F₅ has roots 2 and 3
    Polynomial<F5> p1{F5{73} ,F5{6}, F5{9}, F5{1} };   // for x^2+1
    Polynomial<F5> p2{ F5{1}, F5{1}, F5{1} };   // for x^2+x+1
    std::cout << "Testing f(x)=x^2+1 mod 5:\n";
    for(int i = 0; i < 5; ++i) {
        F5 a{i};
        std::cout << "  f(" << i << ") = " << p1.eval(a) << "\n";
    }
    // You should see f(2)=0 and f(3)=0 → reducible.

    // Now actually divide out one factor:
    p1.reduce();
    std::cout << "After reduce(): degree is now “1” (i.e. linear)\n\n";


    // 2) An irreducible polynomial: x² + x + 1 over F₅ has no roots
    
    std::cout << "Testing g(x)=x^2+x+1 mod 5:\n";
    for(int i = 0; i < 5; ++i) {
        F5 a{i};
        std::cout << "  g(" << i << ") = " << p2.eval(a) << "\n";
    }
    // You should see no zeroes → p2 remains unchanged by reduce().

    return 0;
}
