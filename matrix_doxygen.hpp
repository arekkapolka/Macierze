#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace ak {

/**
 * @brief Base class for matrix operations
 *
 * This template class provides basic matrix operations such as addition, subtraction, multiplication, 
 * scalar multiplication, transposition, and finding the minimum and maximum elements in a matrix.
 * The operations ensure correct validation for matrix dimensions.
 *
 * @tparam T The type of elements stored in the matrix (e.g., double, int, etc.)
 */
template <typename T>
class MatrixBase {
protected:
    size_t m_rows; ///< Number of rows in the matrix
    size_t m_cols; ///< Number of columns in the matrix
    std::vector<std::vector<T>> m_data; ///< The matrix data

public:
    /**
     * @brief Constructor to initialize a matrix with given dimensions.
     * 
     * @param rows The number of rows in the matrix
     * @param cols The number of columns in the matrix
     * 
     * Throws an exception if the dimensions are zero.
     */
    MatrixBase(size_t rows, size_t cols) : m_rows(rows), m_cols(cols), m_data(rows, std::vector<T>(cols)) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Matrix dimensions must be greater than zero");
        }
    }

    /**
     * @brief Accessor to get the element at a specific position.
     * 
     * @param row The row index
     * @param col The column index
     * 
     * @return The reference to the matrix element at the specified position
     * 
     * Throws an exception if the indices are out of range.
     */
    T& at(size_t row, size_t col) {
        if (row >= m_rows || col >= m_cols) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return m_data[row][col];
    }

    /**
     * @brief Const accessor to get the element at a specific position.
     * 
     * @param row The row index
     * @param col The column index
     * 
     * @return The constant reference to the matrix element at the specified position
     * 
     * Throws an exception if the indices are out of range.
     */
    const T& at(size_t row, size_t col) const {
        if (row >= m_rows || col >= m_cols) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return m_data[row][col];
    }

    /**
     * @brief Get the number of rows in the matrix.
     * 
     * @return The number of rows
     */
    size_t getRows() const { return m_rows; }

    /**
     * @brief Get the number of columns in the matrix.
     * 
     * @return The number of columns
     */
    size_t getCols() const { return m_cols; }

    /**
     * @brief Scalar multiplication of the matrix.
     * 
     * @param scalar The scalar to multiply each element of the matrix by
     * @return A reference to the current matrix after scalar multiplication
     */
    MatrixBase& operator*=(T scalar) {
        for (auto& row : m_data) {
            for (auto& element : row) {
                element *= scalar;
            }
        }
        return *this;
    }

    /**
     * @brief Matrix addition.
     * 
     * @param other The matrix to add
     * @return A reference to the current matrix after addition
     * 
     * Throws an exception if the matrices do not have the same dimensions.
     */
    MatrixBase& operator+=(const MatrixBase& other) {
        if (this->m_rows != other.m_rows || this->m_cols != other.m_cols) {
            throw std::invalid_argument("Matrices dimensions must match for addition");
        }

        for (size_t i = 0; i < this->m_rows; ++i) {
            for (size_t j = 0; j < this->m_cols; ++j) {
                this->m_data[i][j] += other.m_data[i][j];
            }
        }
        return *this;
    }

    /**
     * @brief Matrix subtraction.
     * 
     * @param other The matrix to subtract
     * @return A reference to the current matrix after subtraction
     * 
     * Throws an exception if the matrices do not have the same dimensions.
     */
    MatrixBase& operator-=(const MatrixBase& other) {
        if (this->m_rows != other.m_rows || this->m_cols != other.m_cols) {
            throw std::invalid_argument("Matrices dimensions must match for subtraction");
        }

        for (size_t i = 0; i < this->m_rows; ++i) {
            for (size_t j = 0; j < this->m_cols; ++j) {
                this->m_data[i][j] -= other.m_data[i][j];
            }
        }
        return *this;
    }

    /**
     * @brief Matrix multiplication.
     * 
     * @param other The matrix to multiply with
     * @return A new matrix resulting from the multiplication
     * 
     * Throws an exception if the number of columns of the first matrix does not match the number of rows of the second.
     */
    MatrixBase operator*(const MatrixBase& other) const {
        if (this->m_cols != other.m_rows) {
            throw std::invalid_argument("Invalid dimensions for matrix multiplication");
        }

        MatrixBase result(this->m_rows, other.m_cols);
        for (size_t i = 0; i < this->m_rows; ++i) {
            for (size_t j = 0; j < other.m_cols; ++j) {
                result.m_data[i][j] = 0;
                for (size_t k = 0; k < this->m_cols; ++k) {
                    result.m_data[i][j] += this->m_data[i][k] * other.m_data[k][j];
                }
            }
        }
        return result;
    }

    /**
     * @brief Transpose the matrix.
     * 
     * This function swaps the rows and columns of the matrix.
     */
    void transpose() {
        for (size_t i = 0; i < this->m_rows; ++i) {
            for (size_t j = i + 1; j < this->m_cols; ++j) {
                std::swap(this->m_data[i][j], this->m_data[j][i]);
            }
        }
    }

    /**
     * @brief Find the minimum value in the matrix.
     * 
     * @return The minimum value in the matrix
     * 
     * Throws an exception if the matrix is empty.
     */
    T min() const {
        if (this->m_rows == 0 || this->m_cols == 0) {
            throw std::invalid_argument("Cannot find min of an empty matrix");
        }
        T minValue = m_data[0][0];
        for (size_t i = 0; i < m_rows; ++i) {
            for (size_t j = 0; j < m_cols; ++j) {
                if (m_data[i][j] < minValue) {
                    minValue = m_data[i][j];
                }
            }
        }
        return minValue;
    }

    /**
     * @brief Find the maximum value in the matrix.
     * 
     * @return The maximum value in the matrix
     * 
     * Throws an exception if the matrix is empty.
     */
    T max() const {
        if (this->m_rows == 0 || this->m_cols == 0) {
            throw std::invalid_argument("Cannot find max of an empty matrix");
        }
        T maxValue = m_data[0][0];
        for (size_t i = 0; i < m_rows; ++i) {
            for (size_t j = 0; j < m_cols; ++j) {
                if (m_data[i][j] > maxValue) {
                    maxValue = m_data[i][j];
                }
            }
        }
        return maxValue;
    }

    /**
     * @brief Destructor
     * 
     * Virtual destructor to ensure proper cleanup for derived classes.
     */
    virtual ~MatrixBase() = default;
};

/**
 * @brief Derived class for square matrices
 *
 * This class provides additional functionality for square matrices, such as the calculation of determinants.
 *
 * @tparam T The type of elements stored in the square matrix (e.g., double, int, etc.)
 */
template <typename T>
class SquareMatrix : public MatrixBase<T> {
public:
    /**
     * @brief Constructor to initialize a square matrix with the given size.
     * 
     * @param size The size of the square matrix
     * 
     * Throws an exception if the size is zero.
     */
    SquareMatrix(size_t size) : MatrixBase<T>(size, size) {
        if (size == 0) {
            throw std::invalid_argument("Matrix size must be greater than zero");
        }
    }

    /**
     * @brief Calculate the determinant of the square matrix.
     * 
     * This function uses a recursive method to compute the determinant for square matrices.
     * 
     * @return The determinant of the matrix
     */
    T determinant() const {
        if (this->m_rows == 1) {
            return this->m_data[0][0];
        }

        T det = 0;
        for (size_t col = 0; col < this->m_cols; ++col) {
            SquareMatrix subMatrix = this->getSubMatrix(0, col);
            det += (col % 2 == 0 ? 1 : -1) * this->m_data[0][col] * subMatrix.determinant();
        }
        return det;
    }

private:
    /**
     * @brief Helper function to get a submatrix by excluding a specific row and column.
     * 
     * @param excludeRow The row to exclude
     * @param excludeCol The column to exclude
     * 
     * @return A new square matrix that is the submatrix excluding the specified row and column
     */
    SquareMatrix getSubMatrix(size_t excludeRow, size_t excludeCol) const {
        SquareMatrix subMatrix(this->m_rows - 1);
        size_t subRow = 0;
        for (size_t i = 0; i < this->m_rows; ++i) {
            if (i == excludeRow) continue;
            size_t subCol = 0;
            for (size_t j = 0; j < this->m_cols; ++j) {
                if (j == excludeCol) continue;
                subMatrix.at(subRow, subCol++) = this->m_data[i][j];
            }
            ++subRow;
        }
        return subMatrix;
    }
};

/**
 * @brief Overload the output stream operator to print MatrixBase objects.
 * 
 * @param os The output stream
 * @param matrix The matrix object to print
 * @return The output stream with the matrix data
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, const MatrixBase<T>& matrix) {
    for (size_t i = 0; i < matrix.getRows(); ++i) {
        for (size_t j = 0; j < matrix.getCols(); ++j) {
            os << matrix.at(i, j) << " ";
        }
        os << "\n";
    }
    return os;
}

} // namespace ak

#endif // MATRIX_HPP
