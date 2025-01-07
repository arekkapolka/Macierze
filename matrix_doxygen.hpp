#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace ak {

/**
 * @brief A base class for matrix operations.
 * 
 * This class represents a generic matrix with basic operations such as scalar 
 * multiplication, matrix addition, subtraction, multiplication, transposition, 
 * and finding minimum and maximum values.
 * 
 * @tparam T The data type of the matrix elements.
 */
template <typename T>
class MatrixBase {
protected:
    size_t m_rows; ///< The number of rows in the matrix.
    size_t m_cols; ///< The number of columns in the matrix.
    std::vector<std::vector<T>> m_data; ///< The matrix data.

public:
    /**
     * @brief Constructs a matrix with specified dimensions.
     * 
     * Initializes a matrix with the given number of rows and columns.
     * Throws an exception if the dimensions are invalid (zero or negative).
     * 
     * @param rows The number of rows in the matrix.
     * @param cols The number of columns in the matrix.
     * @throws std::invalid_argument If either rows or columns are zero or negative.
     */
    MatrixBase(size_t rows, size_t cols);

    /**
     * @brief Access an element of the matrix.
     * 
     * Provides access to an element at the specified row and column index.
     * 
     * @param row The row index of the element.
     * @param col The column index of the element.
     * @return Reference to the matrix element.
     * @throws std::out_of_range If the row or column index is out of range.
     */
    T& at(size_t row, size_t col);

    /**
     * @brief Access an element of the matrix (constant version).
     * 
     * Provides access to an element at the specified row and column index.
     * 
     * @param row The row index of the element.
     * @param col The column index of the element.
     * @return Constant reference to the matrix element.
     * @throws std::out_of_range If the row or column index is out of range.
     */
    const T& at(size_t row, size_t col) const;

    /**
     * @brief Get the number of rows in the matrix.
     * 
     * @return The number of rows in the matrix.
     */
    size_t getRows() const;

    /**
     * @brief Get the number of columns in the matrix.
     * 
     * @return The number of columns in the matrix.
     */
    size_t getCols() const;

    /**
     * @brief Scalar multiplication of the matrix.
     * 
     * Multiplies each element of the matrix by the given scalar.
     * 
     * @param scalar The scalar value to multiply by.
     * @return Reference to the current matrix after multiplication.
     */
    MatrixBase& operator*=(T scalar);

    /**
     * @brief Matrix addition.
     * 
     * Adds another matrix to the current matrix. The dimensions of both matrices 
     * must match.
     * 
     * @param other The matrix to add.
     * @return Reference to the current matrix after addition.
     * @throws std::invalid_argument If the dimensions of the matrices do not match.
     */
    MatrixBase& operator+=(const MatrixBase& other);

    /**
     * @brief Matrix subtraction.
     * 
     * Subtracts another matrix from the current matrix. The dimensions of both matrices 
     * must match.
     * 
     * @param other The matrix to subtract.
     * @return Reference to the current matrix after subtraction.
     * @throws std::invalid_argument If the dimensions of the matrices do not match.
     */
    MatrixBase& operator-=(const MatrixBase& other);

    /**
     * @brief Matrix multiplication.
     * 
     * Multiplies the current matrix by another matrix. The number of columns in 
     * the first matrix must match the number of rows in the second matrix.
     * 
     * @param other The matrix to multiply by.
     * @return A new matrix containing the result of the multiplication.
     * @throws std::invalid_argument If the matrices' dimensions do not match for multiplication.
     */
    MatrixBase operator*(const MatrixBase& other) const;

    /**
     * @brief Transpose the matrix.
     * 
     * Transposes the current matrix in-place.
     */
    void transpose();

    /**
     * @brief Find the minimum value in the matrix.
     * 
     * Returns the smallest element in the matrix.
     * 
     * @return The minimum value in the matrix.
     * @throws std::invalid_argument If the matrix is empty.
     */
    T min() const;

    /**
     * @brief Find the maximum value in the matrix.
     * 
     * Returns the largest element in the matrix.
     * 
     * @return The maximum value in the matrix.
     * @throws std::invalid_argument If the matrix is empty.
     */
    T max() const;

    /**
     * @brief Virtual destructor.
     */
    virtual ~MatrixBase() = default;
};

/**
 * @brief A derived class for square matrices.
 * 
 * This class represents a square matrix and provides additional functionality 
 * such as calculating the determinant of the matrix.
 * 
 * @tparam T The data type of the matrix elements.
 */
template <typename T>
class SquareMatrix : public MatrixBase<T> {
public:
    /**
     * @brief Constructs a square matrix with the given size.
     * 
     * Initializes a square matrix of the specified size. Throws an exception if 
     * the size is zero.
     * 
     * @param size The size of the square matrix.
     * @throws std::invalid_argument If the size is zero.
     */
    SquareMatrix(size_t size);

    /**
     * @brief Calculate the determinant of the square matrix.
     * 
     * Recursively calculates the determinant of the matrix using cofactor expansion.
     * 
     * @return The determinant of the matrix.
     */
    T determinant() const;

private:
    /**
     * @brief Get a submatrix by excluding a specific row and column.
     * 
     * Generates a submatrix by excluding the given row and column.
     * 
     * @param excludeRow The row to exclude.
     * @param excludeCol The column to exclude.
     * @return A new submatrix.
     */
    SquareMatrix getSubMatrix(size_t excludeRow, size_t excludeCol) const;
};

/**
 * @brief Overload of the stream insertion operator for matrix output.
 * 
 * Outputs the matrix elements to the provided output stream.
 * 
 * @param os The output stream.
 * @param matrix The matrix to output.
 * @return The output stream after inserting the matrix.
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, const MatrixBase<T>& matrix);

} // namespace ak

#endif // MATRIX_HPP
