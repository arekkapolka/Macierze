#include <iostream>
#include <cstdlib>
#include <ctime>
#include "matrix.hpp"
/**
 * @file main.cpp
 * @brief Demonstrates usage of MatrixBase and SquareMatrix classes for matrix operations.
 * 
 * This file demonstrates the functionality of the MatrixBase and SquareMatrix classes 
 * by performing several matrix operations such as addition, subtraction, multiplication, 
 * scalar multiplication, transposition, and calculating the determinant and extreme values 
 * of a square matrix.
 * 
 * The operations are performed on matrices of different sizes (2x2, 5x5), with random numbers 
 * used for the 5x5 matrix.
 * 
 * @note Exceptions are caught and displayed using the standard exception handling mechanism.
 */

int main() {
    try {
        // Create and initialize a 2x2 matrix
        ak::MatrixBase<double> matrix1(2, 2); ///< 2x2 matrix for demonstration.
        matrix1.at(0, 0) = 1.0;
        matrix1.at(0, 1) = 2.0;
        matrix1.at(1, 0) = 3.0;
        matrix1.at(1, 1) = 4.0;

        std::cout << "Matrix 1 (2x2):\n" << matrix1;

        // Create the second 2x2 matrix
        ak::MatrixBase<double> matrix2(2, 2); ///< Second 2x2 matrix for operations.
        matrix2.at(0, 0) = 5.0;
        matrix2.at(0, 1) = 6.0;
        matrix2.at(1, 0) = 7.0;
        matrix2.at(1, 1) = 8.0;

        std::cout << "\nMatrix 2 (2x2):\n" << matrix2;

        // Matrix addition
        matrix1 += matrix2; ///< Perform matrix addition.
        std::cout << "\nMatrix 1 + Matrix 2:\n" << matrix1;

        // Matrix subtraction
        matrix1 -= matrix2; ///< Perform matrix subtraction.
        std::cout << "\nMatrix 1 - Matrix 2:\n" << matrix1;

        // Matrix multiplication
        ak::MatrixBase<double> resultMul = matrix1 * matrix2; ///< Perform matrix multiplication.
        std::cout << "\nMatrix 1 * Matrix 2:\n" << resultMul;

        // Scalar multiplication
        matrix1 *= 2.0; ///< Perform scalar multiplication.
        std::cout << "\nMatrix 1 after scalar multiplication by 2.0:\n" << matrix1;

        // Matrix transposition
        matrix1.transpose(); ///< Perform matrix transposition.
        std::cout << "\nMatrix 1 after transposition:\n" << matrix1;

        // Initialize a 5x5 matrix with random numbers from 1 to 5
        ak::SquareMatrix<int> matrix3(5); ///< 5x5 matrix with random numbers.
        std::srand(static_cast<unsigned int>(std::time(nullptr))); ///< Seed for random number generation.
        for (size_t i = 0; i < matrix3.getRows(); ++i) {
            for (size_t j = 0; j < matrix3.getCols(); ++j) {
                matrix3.at(i, j) = std::rand() % 5 + 1;  ///< Fill matrix with random numbers from 1 to 5.
            }
        }

        std::cout << "\nMatrix 3 (5x5) with random numbers:\n" << matrix3;

        // Calculate the determinant of the 5x5 matrix
        int determinant = matrix3.determinant(); ///< Calculate the determinant of the matrix.
        std::cout << "\nDeterminant of Matrix 3 (5x5): " << determinant << "\n";

        // Find the minimum and maximum values in Matrix 3
        int minVal = matrix3.min(); ///< Find minimum value in the matrix.
        int maxVal = matrix3.max(); ///< Find maximum value in the matrix.
        std::cout << "\nMinimum value in Matrix 3: " << minVal << "\n";
        std::cout << "Maximum value in Matrix 3: " << maxVal << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n"; ///< Error handling for exceptions.
    }

    return 0;
}
