#include <iostream>
#include <vector>
#include <algorithm>    // for std::swap
#include <stdexcept>    // for std::runtime_error

using namespace std;



template<typename T>
vector<vector<T>> matrix_multiply(){



    
}




template<typename T>
vector<T> gauss_jordan(vector<vector<T>>& A,
                       vector<T>& b)
{
    const size_t N = A.size();
    if (b.size() != N)
        throw std::runtime_error("Dimension mismatch");

    // Forward elimination with partial pivoting
    for (size_t i = 0; i < N; ++i) {
        // 1) Read (or re‑read) the pivot
        T pivot = A[i][i];

        // 2) If zero, find a row below with a non‑zero in col i
        if (pivot == T(0)) {
            size_t swapRow = i + 1;
            for (; swapRow < N && A[swapRow][i] == T(0); ++swapRow);
            if (swapRow == N)
                throw std::runtime_error("Matrix is singular");

            std::swap(A[i], A[swapRow]);
            std::swap(b[i], b[swapRow]);
            pivot = A[i][i];  // reload from the swapped‑in row
        }

        // 3) Normalize the pivot row
        for (size_t k = i; k < N; ++k) {
            A[i][k] /= pivot;
        }
        b[i] /= pivot;

        // 4) Eliminate all rows *below* i
        for (size_t j = i + 1; j < N; ++j) {
            T factor = A[j][i];
            for (size_t k = i; k < N; ++k)
                A[j][k] -= factor * A[i][k];
            b[j] -= factor * b[i];
        }
    }

    // Back‑substitution to get the solution x (A is now upper‑triangular with 1s on the diag)
    vector<T> x(N);
    for (int i = int(N) - 1; i >= 0; --i) {
        T sum = T(0);
        for (size_t k = i + 1; k < N; ++k)
            sum += A[i][k] * x[k];
        x[i] = b[i] - sum;
    }
    return x;
}

template<typename T>
void transpose(vector<vector<T>>& A) {
    size_t N = A.size();
    if (N == 0) return;
    size_t M = A[0].size();
    // allocate B as M×N
    vector<vector<T>> B(M, vector<T>(N));
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < M; ++j)
            B[j][i] = A[i][j];
    A = std::move(B); //B is left empty but valid after this transformation
}

template<typename T>
bool element_is_in(const vector<T>& v, const T& elem) {
    for (auto const& x : v)
        if (x == elem) return true;
    return false;
}

template<typename T>
vector<vector<T>> building_matrix_B(const vector<vector<T>>& A, const vector<T>&index) { // This function should return a square matrix
    // Kom ihåg att inte ändra på varken A eller index i det här fallet
    size_t N = A.size();
    size_t M = A[0].size();
    size_t K = index.size();
    vector<vector<T>> B(N, vector<T>(K));

    int k = 0;
    for(size_t j = 0; j < M; j++){
        if (element_is_in(index, j)){
            for (size_t i = 0; i<N; i++){
                B[i][k] = A[i][j]; // fyller upp B matrisen
            }
            k ++;
        }
    }

    return B; //Ger tillbaka det B som behövs för att använda i det här fallet
}


template<typename T>
std::vector<T> pick_elem(
    const std::vector<T>& v,
    const std::vector<size_t>& index,
    bool pick_index = true
) {
    size_t M = v.size();
    // reserve exactly how many we expect
    size_t expect = pick_index
                   ? index.size()
                   : M - index.size();
    std::vector<T> out;
    out.reserve(expect);

    for (size_t i = 0; i < M; ++i) {
        bool in = element_is_in(index, i);
        if ((pick_index && in) || (!pick_index && !in))
            out.push_back(v[i]);
    }
    return out;
}

template<typename T>
std::vector<std::vector<T>>
matrix_multiply(const std::vector<std::vector<T>>& A,
                const std::vector<std::vector<T>>& B)
{
    // 1) Dimension check
    if (A.empty() || B.empty() || A[0].size() != B.size()) {
        throw std::invalid_argument{
            "Incompatible dimensions: A is "
            + std::to_string(A.size()) + "×" + std::to_string(A.empty()?0:A[0].size())
            + ", B is "
            + std::to_string(B.size()) + "×" + std::to_string(B.empty()?0:B[0].size())
        };
    }

    size_t n = A.size();       // # rows of A
    size_t p = A[0].size();    // # cols of A == # rows of B
    size_t m = B[0].size();    // # cols of B

    // 2) Allocate and zero‐initialize C as an n×m matrix:
    //    - Outer vector has n elements.
    //    - Each inner vector has m elements.
    //    - T{} value‐initializes each entry (e.g. 0 for arithmetic types).
    std::vector<std::vector<T>> C(n, std::vector<T>(m, T{}));

    // 3) Standard triple‐loop multiply
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            for (size_t k = 0; k < p; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;  // RVO/move makes this efficient
}





template<typename T>
vector<T> simplex(const vector<vector<T>>&A, const vector<T>&b, const vector<T>&c){ // main simplex method in this case

// I am using the const functionality since i dont want to mutate A, b or c in this case
// what assumtions are we making
// 1. n x m matrix where the extra variables are already implemented in the matrix
// 2. I have to be careful with the overload in this case.
// 3. This also mean i have to tidy upp the code for transportationsnetworks
// 4. The code will return the x values in this case used

const size_t N = A.size();
const size_t M = A[0].size(); // de här delarna får inte ändras
vector<T> index(N); // indexen får inte överstiga antalet delar här

// initilize the starting index
for (size_t i = 0; i < N, i++){
    index[i] = M - N + i; //lagrar de bakersta indexen först

}

// first loop 

vector<vector<T>> B = building_matrix_B(A, index);
auto B_T = B;
vector<T> x = gauss_jordan(B, b); // hämtar ut B delen
transpose(B_T);
vector<T> c_b = pick_elem(c, index);
vector<T> y = gauss_jordan(B_T, c_b);
// i have realised i now actually need a matrix multiplier function in this case 








}




int main() {

return 0;
}
