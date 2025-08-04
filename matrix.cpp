#include <vector>
#include <iostream>


//Roadmao for this class
// It should be able to handle simple matrix calculations like addition, subtraction, multiplication, and division.
// Inverse calculatoin is also needed aswell as // determinant calculation.


class Matrix {
    std::vector<std::vector<double>> data; // 2D vector to hold the matrix data
    public:
    // Constructor to initialize the matrix with a given size and default value
    Matrix(int rows, int cols, double default_value = 0.0) : data(rows, std::vector<double>(cols, default_value)) {}

    // Constructor to initialize the matrix with a 2D vector (MATLAB-like initialization)
    Matrix(const std::vector<std::vector<double>>& input_data) : data(input_data) {}

    // Function to set a value at a specific position in the matrix
};


int main(){

    // Example usage of the Matrix class
    Matrix mat1(3, 3, 1.0); // Create a 3x3 matrix initialized with 1.0
    Matrix mat2({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}); // Create a matrix from a 2D vector

    // Here you can add more operations like addition, subtraction, multiplication etc.

    std::cout << "Matrix operations can be implemented here." << std::endl;

    return 0;       
}