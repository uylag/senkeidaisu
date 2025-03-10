#pragma once
#include <sparse.hpp>
#include <dense.hpp>
#include <variant>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <string>

namespace senkeidaisu 
{
    template <typename TN>
    class Matrix
    {
    public:
        Matrix();

        Matrix(const int64& rows, const int64& cols, 
               const TN& default_value = TN(0));

        Matrix(const TN*& c_matrix, const int64& rows, 
               const int64& cols);

        Matrix(const std::vector<TN>& c_matrix, const int64& rows, 
               const int64& cols);

        Matrix(const std::vector<std::vector<TN>>& matrix);

        Matrix(const Matrix<TN>& other);

        Matrix(Matrix<TN>&& other) noexcept;

        Matrix(std::initializer_list<std::initializer_list<TN>> matrix);

        Matrix(std::initializer_list<TN> matrix, int rows, int columns);

        Matrix(const DenseMatrix<TN>& other);
        Matrix(const SparseMatrix<TN>& other)
        {
            dense_ = DenseMatrix<TN>(other);
            sparse_ = other;
        };

        ~Matrix();

        TN& operator()(int64 i, int64 j);
        const TN& operator()(int64 i, int64 j) const;

        Matrix<TN> operator+(const TN& scalar) const;
        Matrix<TN> operator+(const TN*& c_matrix) const;
        Matrix<TN> operator+(const TN** matrix) const;
        Matrix<TN> operator+(const std::vector<TN> c_matrix) const;
        Matrix<TN> operator+(const std::vector<std::vector<TN>> c_matrix) const;
        Matrix<TN> operator+(const Matrix<TN>& other) const;
        Matrix<TN> operator+(const DenseMatrix<TN>& other) const;
        Matrix<TN> operator+(const SparseMatrix<TN>& other) const;

        Matrix<TN> operator-(const TN& scalar) const;
        Matrix<TN> operator-(const TN*& c_matrix) const;
        Matrix<TN> operator-(const TN** matrix) const;
        Matrix<TN> operator-(const std::vector<TN> c_matrix) const;
        Matrix<TN> operator-(const std::vector<std::vector<TN>> c_matrix) const;
        Matrix<TN> operator-(const Matrix<TN>& other) const;
        Matrix<TN> operator-(const DenseMatrix<TN>& other) const;
        Matrix<TN> operator-(const SparseMatrix<TN>& other) const;

        Matrix<TN> operator*(const TN& scalar) const;
        Matrix<TN> operator*(const TN*& c_matrix) const;
        Matrix<TN> operator*(const TN** matrix) const;
        Matrix<TN> operator*(const std::vector<TN> c_matrix) const;
        Matrix<TN> operator*(const std::vector<std::vector<TN>> c_matrix) const;
        Matrix<TN> operator*(const Matrix<TN>& other) const;
        Matrix<TN> operator*(const DenseMatrix<TN>& other) const;
        Matrix<TN> operator*(const SparseMatrix<TN>& other) const;

        Matrix<TN>& operator=(const Matrix<TN>& other);
        Matrix<TN>& operator=(Matrix<TN>&& other) noexcept;
        Matrix<TN>& operator=(const DenseMatrix<TN>& other) const;
        Matrix<TN>& operator=(const SparseMatrix<TN>& other) const;

        Matrix<TN> operator-() const;

        Matrix<TN>& operator+=(const TN& scalar);
        Matrix<TN>& operator+=(const TN*& c_matrix);
        Matrix<TN>& operator+=(const TN** matrix);
        Matrix<TN>& operator+=(const std::vector<TN> c_matrix);
        Matrix<TN>& operator+=(const std::vector<std::vector<TN>> c_matrix);
        Matrix<TN>& operator+=(const Matrix<TN>& other);
        Matrix<TN>& operator+=(const DenseMatrix<TN>& other) const;
        Matrix<TN>& operator+=(const SparseMatrix<TN>& other) const;

        Matrix<TN>& operator-=(const TN& scalar);
        Matrix<TN>& operator-=(const TN*& c_matrix);
        Matrix<TN>& operator-=(const TN** matrix);
        Matrix<TN>& operator-=(const std::vector<TN> c_matrix);
        Matrix<TN>& operator-=(const std::vector<std::vector<TN>> c_matrix);
        Matrix<TN>& operator-=(const Matrix<TN>& other);
        Matrix<TN>& operator-=(const DenseMatrix<TN>& other) const;
        Matrix<TN>& operator-=(const SparseMatrix<TN>& other) const;

        Matrix<TN>& operator*=(const TN& scalar);
        Matrix<TN>& operator*=(const TN*& c_matrix);
        Matrix<TN>& operator*=(const TN** matrix);
        Matrix<TN>& operator*=(const std::vector<TN> c_matrix);
        Matrix<TN>& operator*=(const std::vector<std::vector<TN>> c_matrix);
        Matrix<TN>& operator*=(const Matrix<TN>& other);
        Matrix<TN>& operator*=(const DenseMatrix<TN>& other) const;
        Matrix<TN>& operator*=(const SparseMatrix<TN>& other) const;

        Matrix<TN> operator/(const TN& scalar) const;
        Matrix<TN>& operator/=(const TN& scalar);

        bool operator==(const Matrix<TN>& other) const { return this->dense_ == other.dense_; }
        bool operator!=(const Matrix<TN>& other) const { return !(*this == other); }

        friend std::ostream& operator<<(std::ostream& os, const Matrix<TN>& matrix)
        {
            for (int64 i = 0; i < matrix.rows(); ++i)
            {
                for (int64 j = 0; j < matrix.cols(); ++j)
                {
                    os << matrix(i, j) << " ";
                }
                os << '\n';
            }
            return os;
        }

        inline std::vector<std::vector<TN>> dd() const { return dense_.dd(); }
        inline std::vector<TN>              v() const { return dense_.v(); }

        inline Matrix<TN> T(bool inplace = false);

        inline bool is_empty() const;
        inline bool is_square() const;

        inline int64 rows() const;
        inline static int64 rows(const Matrix& matrix) { return matrix.rows(); }

        inline int64 cols() const;
        inline static int64 cols(const Matrix& matrix) { return matrix.cols(); }

        inline int64 size() const;
        inline static int64 size(const Matrix& matrix) { return matrix.size(); }

        inline static Matrix<TN> identity(int64 size);
        inline static Matrix<TN> diagonal(const std::vector<TN>& values);

        inline Matrix<TN> submatrix(const std::initializer_list<int>& rows_start_end,
                                    const std::initializer_list<int>& columns_start_end);

        inline void fill(const TN& value, bool remain_filled = true);
        inline void resize(int64 rows, int64 cols);

        Matrix<TN> clone() const;

        static Matrix<TN> read_csv(const std::string& location)
        {
            // Stub
            return Matrix<TN>();
        };

    private:
        SparseMatrix<TN> sparse_;
        DenseMatrix<TN> dense_;
    };

    /////////////// IMPLEMENTATION ///////////////////

    // Constructors

    template <typename TN>
    Matrix<TN>::Matrix() : dense_(), sparse_() {}

    template <typename TN>
    Matrix<TN>::Matrix(const int64 &rows, const int64 &cols, const TN &default_value)
        : dense_(rows, cols, default_value), sparse_(rows, cols, default_value) {}

    template <typename TN>
    Matrix<TN>::Matrix(const TN*& c_matrix, const int64 &rows, const int64 &cols)
        : dense_(c_matrix, rows, cols), sparse_(c_matrix, rows, cols) {}

    template <typename TN>
    Matrix<TN>::Matrix(const std::vector<TN>& c_matrix, const int64 &rows, const int64 &cols)
        : dense_(c_matrix, rows, cols), sparse_(c_matrix, rows, cols) {}

    template <typename TN>
    Matrix<TN>::Matrix(const std::vector<std::vector<TN>>& matrix)
        : dense_(matrix), sparse_(matrix) {}

    template <typename TN>
    Matrix<TN>::Matrix(const Matrix<TN>& other)
        : dense_(other.dense_), sparse_(other.sparse_) {}

    template <typename TN>
    Matrix<TN>::Matrix(Matrix<TN>&& other) noexcept
        : dense_(std::move(other.dense_)), sparse_(std::move(other.sparse_)) {}

    template <typename TN>
    Matrix<TN>::Matrix(std::initializer_list<std::initializer_list<TN>> matrix)
        : dense_(matrix), sparse_(matrix) {}

    template <typename TN>
    Matrix<TN>::Matrix(std::initializer_list<TN> matrix, int rows, int columns)
        : dense_(matrix, rows, columns),
          sparse_(std::vector<TN>(matrix), rows, columns) {}

    template <typename TN>
    Matrix<TN>::Matrix(const DenseMatrix<TN>& other)
        : dense_(other), sparse_(other) {}

    template <typename TN>
    Matrix<TN>::~Matrix() {}

    // Element Access

    template <typename TN>
    TN& Matrix<TN>::operator()(int64 i, int64 j) {
        return dense_(i, j);
    }

    template <typename TN>
    const TN& Matrix<TN>::operator()(int64 i, int64 j) const {
        return dense_(i, j);
    }

    // Assignment Operators

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator=(const Matrix<TN>& other) {
        if (this == &other)
            return *this;
        dense_ = other.dense_;
        sparse_ = other.sparse_;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator=(Matrix<TN>&& other) noexcept {
        if (this == &other)
            return *this;
        dense_ = std::move(other.dense_);
        sparse_ = std::move(other.sparse_);
        return *this;
    }

    // Assignment from DenseMatrix and SparseMatrix (declared const; using const_cast)
    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator=(const DenseMatrix<TN>& other) const {
        const_cast<DenseMatrix<TN>&>(dense_) = other;
        const_cast<SparseMatrix<TN>&>(sparse_) = SparseMatrix<TN>(other);
        return const_cast<Matrix<TN>&>(*this);
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator=(const SparseMatrix<TN>& other) const {
        const_cast<DenseMatrix<TN>&>(dense_) = DenseMatrix<TN>(other);
        const_cast<SparseMatrix<TN>&>(sparse_) = other;
        return const_cast<Matrix<TN>&>(*this);
    }

    // Arithmetic Operators (for scalar and Matrix)

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator+(const TN& scalar) const {
        DenseMatrix<TN> result_dense = dense_ + scalar;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator+(const TN*& c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ + c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator+(const TN** matrix) const {
        DenseMatrix<TN> result_dense = dense_ + matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator+(const std::vector<TN> c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ + c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator+(const std::vector<std::vector<TN>> c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ + c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator+(const Matrix<TN>& other) const {
        if (dense_.rows() != other.dense_.rows() || dense_.cols() != other.dense_.cols())
            throw std::invalid_argument("Matrix dimensions do not match for addition");
        DenseMatrix<TN> result_dense = dense_ + other.dense_;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator+(const DenseMatrix<TN>& other) const {
        DenseMatrix<TN> result_dense = dense_ + other;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator+(const SparseMatrix<TN>& other) const {
        DenseMatrix<TN> temp(other);
        DenseMatrix<TN> result_dense = dense_ + temp;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-(const TN& scalar) const {
        DenseMatrix<TN> result_dense = dense_ - scalar;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-(const TN*& c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ - c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-(const TN** matrix) const {
        DenseMatrix<TN> result_dense = dense_ - matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-(const std::vector<TN> c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ - c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-(const std::vector<std::vector<TN>> c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ - c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-(const Matrix<TN>& other) const {
        if (dense_.rows() != other.dense_.rows() || dense_.cols() != other.dense_.cols())
            throw std::invalid_argument("Matrix dimensions do not match for subtraction");
        DenseMatrix<TN> result_dense = dense_ - other.dense_;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-(const DenseMatrix<TN>& other) const {
        DenseMatrix<TN> result_dense = dense_ - other;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-(const SparseMatrix<TN>& other) const {
        DenseMatrix<TN> temp(other);
        DenseMatrix<TN> result_dense = dense_ - temp;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator*(const TN& scalar) const {
        DenseMatrix<TN> result_dense = dense_ * scalar;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator*(const TN*& c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ * c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator*(const TN** matrix) const {
        DenseMatrix<TN> result_dense = dense_ * matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator*(const std::vector<TN> c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ * c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator*(const std::vector<std::vector<TN>> c_matrix) const {
        DenseMatrix<TN> result_dense = dense_ * c_matrix;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator*(const Matrix<TN>& other) const {
        if (dense_.cols() != other.dense_.rows())
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        SparseMatrix<TN> result_sparse = sparse_ * other.sparse_;
        DenseMatrix<TN> result_dense(result_sparse);
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator*(const DenseMatrix<TN>& other) const {
        if (dense_.cols() != other.rows())
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        DenseMatrix<TN> result_dense = dense_ * other;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator*(const SparseMatrix<TN>& other) const {
        SparseMatrix<TN> result_sparse = sparse_ * other;
        return Matrix<TN>(result_sparse);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator-() const {
        DenseMatrix<TN> result_dense = -dense_;
        return Matrix<TN>(result_dense);
    }

    // Compound Assignment Operators for scalar and Matrix

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator+=(const TN& scalar) {
        dense_ += scalar;
        sparse_ = SparseMatrix<TN>(dense_);
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator+=(const TN*& c_matrix) {
        *this = *this + c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator+=(const TN** matrix) {
        *this = *this + matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator+=(const std::vector<TN> c_matrix) {
        *this = *this + c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator+=(const std::vector<std::vector<TN>> c_matrix) {
        *this = *this + c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator+=(const Matrix<TN>& other) {
        if (dense_.rows() != other.dense_.rows() || dense_.cols() != other.dense_.cols())
            throw std::invalid_argument("Matrix dimensions do not match for addition");
        dense_ += other.dense_;
        sparse_ = SparseMatrix<TN>(dense_);
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator+=(const DenseMatrix<TN>& other) const {
        const_cast<DenseMatrix<TN>&>(dense_) += other;
        const_cast<SparseMatrix<TN>&>(sparse_) = SparseMatrix<TN>(dense_);
        return const_cast<Matrix<TN>&>(*this);
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator+=(const SparseMatrix<TN>& other) const {
        DenseMatrix<TN> temp(other);
        const_cast<DenseMatrix<TN>&>(dense_) += temp;
        const_cast<SparseMatrix<TN>&>(sparse_) = SparseMatrix<TN>(dense_);
        return const_cast<Matrix<TN>&>(*this);
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator-=(const TN& scalar) {
        dense_ -= scalar;
        sparse_ = SparseMatrix<TN>(dense_);
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator-=(const TN*& c_matrix) {
        *this = *this - c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator-=(const TN** matrix) {
        *this = *this - matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator-=(const std::vector<TN> c_matrix) {
        *this = *this - c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator-=(const std::vector<std::vector<TN>> c_matrix) {
        *this = *this - c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator-=(const Matrix<TN>& other) {
        if (dense_.rows() != other.dense_.rows() || dense_.cols() != other.dense_.cols())
            throw std::invalid_argument("Matrix dimensions do not match for subtraction");
        dense_ -= other.dense_;
        sparse_ = SparseMatrix<TN>(dense_);
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator-=(const DenseMatrix<TN>& other) const {
        const_cast<DenseMatrix<TN>&>(dense_) -= other;
        const_cast<SparseMatrix<TN>&>(sparse_) = SparseMatrix<TN>(dense_);
        return const_cast<Matrix<TN>&>(*this);
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator-=(const SparseMatrix<TN>& other) const {
        DenseMatrix<TN> temp(other);
        const_cast<DenseMatrix<TN>&>(dense_) -= temp;
        const_cast<SparseMatrix<TN>&>(sparse_) = SparseMatrix<TN>(dense_);
        return const_cast<Matrix<TN>&>(*this);
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator*=(const TN& scalar) {
        dense_ *= scalar;
        sparse_ = SparseMatrix<TN>(dense_);
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator*=(const TN*& c_matrix) {
        *this = *this * c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator*=(const TN** matrix) {
        *this = *this * matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator*=(const std::vector<TN> c_matrix) {
        *this = *this * c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator*=(const std::vector<std::vector<TN>> c_matrix) {
        *this = *this * c_matrix;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator*=(const Matrix<TN>& other) {
        *this = *this * other;
        return *this;
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator*=(const DenseMatrix<TN>& other) const {
        *this = *this * other;
        return const_cast<Matrix<TN>&>(*this);
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator*=(const SparseMatrix<TN>& other) const {
        *this = *this * other;
        return const_cast<Matrix<TN>&>(*this);
    }

    // Division Operators

    template <typename TN>
    Matrix<TN> Matrix<TN>::operator/(const TN& scalar) const {
        if (scalar == TN(0))
            throw std::invalid_argument("Division by zero");
        DenseMatrix<TN> result_dense = dense_ / scalar;
        return Matrix<TN>(result_dense);
    }

    template <typename TN>
    Matrix<TN>& Matrix<TN>::operator/=(const TN& scalar) {
        if (scalar == TN(0))
            throw std::invalid_argument("Division by zero");
        dense_ /= scalar;
        sparse_ = SparseMatrix<TN>(dense_);
        return *this;
    }

    // Utility Functions

    template <typename TN>
    Matrix<TN> Matrix<TN>::T(bool inplace) {
        DenseMatrix<TN> transposed = dense_.T();
        if (inplace) {
            dense_ = std::move(transposed);
            sparse_ = SparseMatrix<TN>(dense_);
            return *this;
        }
        return Matrix<TN>(transposed);
    }

    template <typename TN>
    bool Matrix<TN>::is_empty() const {
        return size() == 0;
    }

    template <typename TN>
    bool Matrix<TN>::is_square() const {
        return rows() == cols();
    }

    template <typename TN>
    int64 Matrix<TN>::rows() const {
        return dense_.rows();
    }

    template <typename TN>
    int64 Matrix<TN>::cols() const {
        return dense_.cols();
    }

    template <typename TN>
    int64 Matrix<TN>::size() const {
        return rows() * cols();
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::identity(int64 size) {
        DenseMatrix<TN> id = DenseMatrix<TN>::identity(size);
        return Matrix<TN>(id);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::diagonal(const std::vector<TN>& values) {
        DenseMatrix<TN> diag = DenseMatrix<TN>::diagonal(values);
        return Matrix<TN>(diag);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::submatrix(const std::initializer_list<int>& rows_range,
                                     const std::initializer_list<int>& cols_range) {
        if (rows_range.size() != 2 || cols_range.size() != 2)
            throw std::invalid_argument("Submatrix ranges must have exactly 2 elements (start and end)");
        auto rs = *rows_range.begin();
        auto re = *(rows_range.begin() + 1);
        auto cs = *cols_range.begin();
        auto ce = *(cols_range.begin() + 1);
        DenseMatrix<TN> sub = dense_.submatrix({rs, re}, {cs, ce});
        return Matrix<TN>(sub);
    }

    template <typename TN>
    void Matrix<TN>::fill(const TN& value, bool remain_filled) {
        dense_.fill(value, remain_filled);
        sparse_ = SparseMatrix<TN>(dense_);
    }

    template <typename TN>
    void Matrix<TN>::resize(int64 rows, int64 cols) {
        dense_.resize(rows, cols);
        sparse_ = SparseMatrix<TN>(dense_);
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::clone() const {
        return Matrix<TN>(dense_);
    }
}