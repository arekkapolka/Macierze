#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace ak {

// Base class for matrix operations
template <typename T>
class MatrixBase {
protected:
    size_t m_rows; // Number of rows
    size_t m_cols; // Number of columns
    std::vector<std::vector<T>> m_data; // Matrix data

public:
    MatrixBase(size_t rows, size_t cols) : m_rows(rows), m_cols(cols), m_data(rows, std::vector<T>(cols)) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Matrix dimensions must be greater than zero");
        }
    }

    T& at(size_t row, size_t col) {
        if (row >= m_rows || col >= m_cols) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return m_data[row][col];
    }

    const T& at(size_t row, size_t col) const {
        if (row >= m_rows || col >= m_cols) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return m_data[row][col];
    }

    size_t getRows() const { return m_rows; }
    size_t getCols() const { return m_cols; }

    // Scalar multiplication
    MatrixBase& operator*=(T scalar) {
        for (auto& row : m_data) {
            for (auto& element : row) {
                element *= scalar;
            }
        }
        return *this;
    }

    // Matrix addition
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

    // Matrix subtraction
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

    // Matrix multiplication
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

    // Transpose of the matrix
    void transpose() {
        for (size_t i = 0; i < this->m_rows; ++i) {
            for (size_t j = i + 1; j < this->m_cols; ++j) {
                std::swap(this->m_data[i][j], this->m_data[j][i]);
            }
        }
    }

    // Find minimum value in the matrix
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

    // Find maximum value in the matrix
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

    virtual ~MatrixBase() = default;
};

// Derived class for square matrices
template <typename T>
class SquareMatrix : public MatrixBase<T> {
public:
    SquareMatrix(size_t size) : MatrixBase<T>(size, size) {
        if (size == 0) {
            throw std::invalid_argument("Matrix size must be greater than zero");
        }
    }

    // Calculate determinant (recursive method)
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

// Overload the operator<< to print MatrixBase objects
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
