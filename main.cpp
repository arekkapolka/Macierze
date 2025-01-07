#include <iostream>
#include <cstdlib>
#include <ctime>
#include "matrix.hpp"

int main() {
    try {
        // Create and initialize a 2x2 matrix
        ak::MatrixBase<double> matrix1(2, 2);
        matrix1.at(0, 0) = 1.0;
        matrix1.at(0, 1) = 2.0;
        matrix1.at(1, 0) = 3.0;
        matrix1.at(1, 1) = 4.0;

        std::cout << "Matrix 1 (2x2):\n" << matrix1;

        // Create the second 2x2 matrix
        ak::MatrixBase<double> matrix2(2, 2);
        matrix2.at(0, 0) = 5.0;
        matrix2.at(0, 1) = 6.0;
        matrix2.at(1, 0) = 7.0;
        matrix2.at(1, 1) = 8.0;

        std::cout << "\nMatrix 2 (2x2):\n" << matrix2;

        // Matrix addition
        matrix1 += matrix2;
        std::cout << "\nMatrix 1 + Matrix 2:\n" << matrix1;

        // Matrix subtraction
        matrix1 -= matrix2;
        std::cout << "\nMatrix 1 - Matrix 2:\n" << matrix1;

        // Matrix multiplication
        ak::MatrixBase<double> resultMul = matrix1 * matrix2;
        std::cout << "\nMatrix 1 * Matrix 2:\n" << resultMul;

        // Scalar multiplication
        matrix1 *= 2.0;
        std::cout << "\nMatrix 1 after scalar multiplication by 2.0:\n" << matrix1;

        // Matrix transposition
        matrix1.transpose();
        std::cout << "\nMatrix 1 after transposition:\n" << matrix1;

        // Initialize a 5x5 matrix with random numbers from 1 to 5
        ak::SquareMatrix<int> matrix3(5);
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        for (size_t i = 0; i < matrix3.getRows(); ++i) {
            for (size_t j = 0; j < matrix3.getCols(); ++j) {
                matrix3.at(i, j) = std::rand() % 5 + 1;  // Random number from 1 to 5
            }
        }

        std::cout << "\nMatrix 3 (5x5) with random numbers:\n" << matrix3;

        // Calculate the determinant of the 5x5 matrix
        int determinant = matrix3.determinant();
        std::cout << "\nDeterminant of Matrix 3 (5x5): " << determinant << "\n";

        // Find the minimum and maximum values in Matrix 3
        int minVal = matrix3.min();
        int maxVal = matrix3.max();
        std::cout << "\nMinimum value in Matrix 3: " << minVal << "\n";
        std::cout << "Maximum value in Matrix 3: " << maxVal << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return 0;
}
