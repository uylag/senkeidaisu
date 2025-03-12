#pragma once
#include <sparse.hpp>
#include <dense.hpp>
#include <variant>
#include <cmath>
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

        inline Matrix<TN> submatrix(const std::initializer_list<int64>& rows_start_end,
                                    const std::initializer_list<int64>& columns_start_end);

        inline void       fill(const TN& value, bool remain_filled = true);
        inline void       resize(int64 rows, int64 cols);

        Matrix<TN>        clone() const;

        static Matrix<TN> read_csv(const std::string& location)
        {
            // Stub
            return Matrix<TN>();
        };

        Matrix<TN>        swap_rows(const int64& idx1, const int64& idx2, bool inplace = true);

        Matrix<TN>        ref(bool inplace = true);
        Matrix<TN>        rref(bool inplace = true);

        const TN          det() const;

        Matrix<TN>        LU() const;
        Matrix<TN>        get_L();
        Matrix<TN>        get_U();

        std::pair<Matrix<TN>, std::vector<TN>>        Householder_QR(const Matrix<TN>& A);
        Matrix<TN>                                    QR();
        Matrix<TN>                                    get_R();
        Matrix<TN>                                    get_Q();

        Matrix<TN>                                    get_col(const int64& idx)
        {
            int64 num_rows = this->rows();
            Matrix<TN> result(num_rows, 1, 0);
            for (int64 i = 0; i < num_rows; ++i)
            {
                result(i, idx) = (*this)(i, idx);
            }

            return result;
        };

        int64              sign(TN value);
        long double        norm(int64 p = 2)
        {
            long double sum = 0;
            for (int64 i = 0; i < this->rows(); ++i)
            {
                for (int64 j = 0; j < this->cols(); ++j)
                {
                    sum += std::pow(std::abs((*this)(i, j)), p);
                }
            }
            return std::pow(sum, 1.0 / p);
        }

        Matrix<TN> inv();

    private:
        SparseMatrix<TN> sparse_;
        DenseMatrix<TN> dense_;
        DenseMatrix<TN> QH_;
        std::vector<TN> betas_;
        std::vector<Matrix<TN>> vs_;
    };

    /////////////// IMPLEMENTATION ///////////////////

    // # Math: A^{-1} = R^{-1}Q^T
    // # Math: Q^{-1} = Q^T
    template <typename TN>
    Matrix<TN> Matrix<TN>::inv() {
        if (rows() != cols()) {
            throw std::invalid_argument("Matrix must be square to compute inverse.");
        }

        Matrix<TN> R = this->QR();
        Matrix<TN> Q = this->get_Q();

        // Check for singularity
        int64 n = R.rows();
        const TN tolerance = TN(1e-10);
        for (int64 i = 0; i < n; ++i) {
            if (std::abs(R(i, i)) < tolerance) {
                throw std::runtime_error("Matrix is singular or nearly singular and cannot be inverted.");
            }
        }

        Matrix<TN> Q_T = Q.T(false);
        Matrix<TN> R_inv(n, n, TN(0));

        for (int64 col = 0; col < n; ++col) {
            Matrix<TN> e(n, 1, TN(0));
            e(col, 0) = TN(1);
            for (int64 i = n - 1; i >= 0; --i) {
                TN sum = TN(0);
                for (int64 j = i + 1; j < n; ++j) {
                    sum += R(i, j) * R_inv(j, col);
                }
                R_inv(i, col) = (e(i, 0) - sum) / R(i, i);
            }
        }

        return R_inv * Q_T;
    }


    template <typename TN>
    int64 Matrix<TN>::sign(TN value)
    {
        return (TN(0) < value) - (value < TN(0));
    }

    // # Math: H = I - 2vv^T
    // # Math: v = a - \alpha e_1
    // # Math: e_1 \text{ is the first standard basis vector } (1, 0, \dots, 0)^T 
    template <typename TN>
    std::pair<Matrix<TN>, std::vector<TN>> Matrix<TN>::Householder_QR(const Matrix<TN>& A) {
        int64 num_rows = A.rows();
        int64 num_cols = A.cols();
        int64 k = std::min(num_rows, num_cols);
        Matrix<TN> QH(A);  // Work on a copy of A
        std::vector<TN> betas(k, TN(0));
        vs_.clear();  // Clear previous vectors
        const TN epsilon = 1e-10;

        auto sign_lambda = [](TN x) -> TN {
            return (x >= TN(0)) ? TN(1) : TN(-1);
        };

        for (int64 j = 0; j < k; ++j) {
            Matrix<TN> x = QH.submatrix({j, num_rows}, {j, j + 1});
            TN norm_x = x.norm();
            if (std::abs(norm_x) < epsilon) {
                betas[j] = TN(0);
                continue;
            }
            TN alpha = -sign_lambda(QH(j, j)) * norm_x;
            Matrix<TN> v(x);
            v(0, 0) -= alpha;
            TN beta = (v.T(false) * v)(0, 0);
            if (std::abs(beta) < epsilon) {
                betas[j] = TN(0);
                continue;
            }
            betas[j] = TN(2) / beta;
            vs_.push_back(v);  // Store the Householder vector

            for (int64 col = j; col < num_cols; ++col) {
                Matrix<TN> subCol = QH.submatrix({j, num_rows}, {col, col + 1});
                TN dot = (v.T(false) * subCol)(0, 0);
                TN factor = dot * betas[j];
                Matrix<TN> update = v * factor;
                for (int64 row = j; row < num_rows; ++row) {
                    QH(row, col) -= update(row - j, 0);
                }
            }
        }
        return {QH, betas};
    }

    // QR(): Return the upper triangular R factor.
    template <typename TN>
    Matrix<TN> Matrix<TN>::QR() {
        auto qr_result = Householder_QR(*this);
        QH_ = qr_result.first.dense_;  // Store QH
        betas_ = qr_result.second;  // Store betas
        int64 num_rows = this->rows();
        int64 num_cols = this->cols();
        Matrix<TN> R(num_rows, num_cols, TN(0));
        for (int64 i = 0; i < num_rows; ++i) {
            for (int64 j = i; j < num_cols; ++j) {
                R(i, j) = QH_(i, j);
            }
        }
        return R;
    }

    // get_R(): Same as QR() – returns the R factor.
    template <typename TN>
    Matrix<TN> Matrix<TN>::get_R() {
        return this->QR();
    }

    // get_Q(): Reconstruct the orthogonal Q matrix from the Householder vectors.
    template <typename TN>
    Matrix<TN> Matrix<TN>::get_Q() 
    {
        int64 num_rows = this->rows();
        int64 num_cols = this->cols();
        int64 k = std::min(num_rows, num_cols);
        Matrix<TN> Q = Matrix<TN>::identity(num_rows);
        std::cout << "Identity Q dimensions: " << Q.rows() << "x" << Q.cols() << std::endl;

        // Check if Householder vectors were computed
        if (vs_.size() < static_cast<size_t>(k)) {
            std::cerr << "Warning: vs_ vector size (" << vs_.size() 
                    << ") is less than expected (" << k << ")." << std::endl;
        }

        for (int64 j = k - 1; j >= 0; --j) {
            if (j >= static_cast<int64>(vs_.size())) {
                std::cerr << "Error: vs_ vector is missing element at index " << j << std::endl;
                continue;
            }
            Matrix<TN> v = vs_[j];  // Use stored Householder vector
            // Debug: print v dimensions
            std::cout << "vs_[" << j << "] dimensions: " << v.rows() << "x" << v.cols() << std::endl;
            Matrix<TN> Q_sub = Q.submatrix({ j, Q.rows() }, { 0, Q.cols() });
            std::cout << "Submatrix Q from row " << j << " to " << Q.rows() 
                    << ", all columns, dimensions: " << Q_sub.rows() << "x" << Q_sub.cols() << std::endl;
            Matrix<TN> vT = v.T(false);
            std::cout << "vT dimensions: " << vT.rows() << "x" << vT.cols() << std::endl;
            Matrix<TN> prod = vT * Q_sub;
            std::cout << "Product (vT*Q_sub) dimensions: " << prod.rows() << "x" << prod.cols() << std::endl;
            Matrix<TN> temp = prod * betas_[j];
            std::cout << "After multiplying with beta (" << betas_[j] << "), temp dimensions: " 
                    << temp.rows() << "x" << temp.cols() << std::endl;
            Matrix<TN> update = v * temp;
            std::cout << "Update dimensions: " << update.rows() << "x" << update.cols() << std::endl;
            for (int64 i = j; i < num_rows; ++i) {
                for (int64 col = 0; col < Q.cols(); ++col) {
                    Q(i, col) -= update(i - j, col);
                }
            }
        }
        return Q;
    }


    template <typename TN>
    const TN Matrix<TN>::det() const
    {
        if (!is_square()) {
            throw std::invalid_argument("Determinant is defined only for square matrices.");
        }
        Matrix<TN> lu = this->LU();
        TN d = TN(1);
        int64 n = lu.rows();
        for (int64 i = 0; i < n; ++i) {
            d *= lu(i, i);
        }
        return d;
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::get_L()
    {
        int64 num_rows = this->rows();
        int64 num_cols = this->cols();

        if (num_rows != num_cols)
            throw std::invalid_argument("L can only be extracted from a square LU matrix!");

        Matrix<TN> L(num_rows, num_cols, TN(0));  // Initialize with zeros

        for (int64 i = 0; i < num_rows; ++i)
        {
            L(i, i) = TN(1);  // Diagonal elements of L are 1
            for (int64 j = 0; j < i; ++j)
            {
                L(i, j) = (*this)(i, j);  // Copy lower triangular elements from LU matrix
            }
        }

        return L;
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::get_U()
    {
        int64 num_rows = this->rows();
        int64 num_cols = this->cols();

        if (num_rows != num_cols)
            throw std::invalid_argument("U can only be extracted from a square LU matrix!");

        Matrix<TN> U(num_rows, num_cols, TN(0));  // Initialize with zeros

        for (int64 i = 0; i < num_rows; ++i)
        {
            for (int64 j = i; j < num_cols; ++j)
            {
                U(i, j) = (*this)(i, j);  // Copy upper triangular elements from LU matrix
            }
        }

        return U;
    }

    // # Math: A = LU
    // # Math: L = \begin{bmatrix} 1 & 0 & 0 & \cdots & 0 \\ l_{21} & 1 & 0 & \cdots & 0 \\ l_{31} & l_{32} & 1 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ l_{n1} & l_{n2} & l_{n3} \cdots & 1 \end{bmatrix}
    // # Math: u_{ij} = A_{ij} - \sum_{k=1}^{j-1}l_{ik}u_{kj} \forall i (i \le j)
    // # Math: l_{ij} = \frac{A_{ij} - \sum_{k=1}^{j-1}l_{ik}u_{kj}}{u_{ii}} \forall i (i > j)
    template <typename TN>
    Matrix<TN> Matrix<TN>::LU() const {
        int64 num_rows = this->rows();
        int64 num_cols = this->cols();
        if (num_rows != num_cols)
            throw std::invalid_argument("LU Decomposition requires a square matrix");

        Matrix<TN> result(*this);
        for (int i = 0; i < num_rows; ++i) {
            for (int j = i; j < num_cols; ++j) {
                TN sum = 0;
                for (int k = 0; k < i; ++k) {
                    sum += result(i, k) * result(k, j);
                }
                result(i, j) -= sum;
            }
            for (int j = i + 1; j < num_rows; ++j) {
                TN sum = 0;
                for (int k = 0; k < i; ++k) {
                    sum += result(j, k) * result(k, i);
                }
                if (result(i, i) == TN(0))
                    throw std::runtime_error("Matrix is singular and cannot be decomposed");
                result(j, i) = (result(j, i) - sum) / result(i, i);
            }
        }
        return result;
    }


    template <typename TN>
    Matrix<TN> Matrix<TN>::swap_rows(const int64& idx1, const int64& idx2, bool inplace)
    {
        int64 ind1 = idx1, ind2 = idx2;

        if (idx1 < 0) ind1 += rows();
        if (idx2 < 0) ind2 += rows();

        if (ind1 < 0 || ind1 >= rows() || ind2 < 0 || ind2 >= rows())
            throw std::out_of_range("Row indices out of bounds");

        Matrix<TN>* target = this;
        Matrix<TN> result;

        if (!inplace)
        {
            result = *this;
            target = &result;
        }

        for (int64 j = 0; j < target->cols(); ++j)
            std::swap((*target)(ind1, j), (*target)(ind2, j));

        return *target;
    }

    template <typename TN>
    Matrix<TN> Matrix<TN>::rref(bool inplace)
    {
        // If the matrix is empty, just return itself.
        if (this->rows() == 0) return *this;

        Matrix<TN>* target = this;
        Matrix<TN> result;

        // If not modifying in place, work on a copy.
        if (!inplace)
        {
            result = *this;
            target = &result;
        }

        int64 pivot_row = 0;
        int64 num_rows = target->rows();
        int64 num_cols = target->cols();

        // Process each column.
        for (int64 col = 0; col < num_cols && pivot_row < num_rows; ++col)
        {
            // Find the pivot element in the current column.
            int64 pivot_idx = -1;
            TN max_val = TN(0);
            for (int64 row = pivot_row; row < num_rows; ++row)
            {
                if (std::abs((*target)(row, col)) > std::abs(max_val))
                {
                    max_val = (*target)(row, col);
                    pivot_idx = row;
                }
            }

            // If no non-zero pivot found in this column, move to next column.
            if (pivot_idx == -1 || max_val == TN(0))
                continue;

            // Swap the current pivot row with the row having the maximum pivot if needed.
            if (pivot_idx != pivot_row)
            {
                target->swap_rows(pivot_row, pivot_idx, true);
            }

            // Normalize the pivot row by dividing by the pivot element.
            TN pivot_value = (*target)(pivot_row, col);
            for (int64 j = col; j < num_cols; ++j)
            {
                (*target)(pivot_row, j) /= pivot_value;
            }

            // Eliminate all other entries in the pivot column.
            for (int64 row = 0; row < num_rows; ++row)
            {
                if (row != pivot_row)
                {
                    TN factor = (*target)(row, col);
                    for (int64 j = col; j < num_cols; ++j)
                    {
                        (*target)(row, j) -= factor * (*target)(pivot_row, j);
                    }
                }
            }

            // Move to the next pivot row.
            pivot_row++;
        }

        return *target;
    }


    template <typename TN>
    Matrix<TN> Matrix<TN>::ref(bool inplace) 
    {
        if (this->rows() == 1 || this->rows() == 0) return *this;

        Matrix<TN>* target = this;
        Matrix<TN> result;

        if (!inplace) 
        {
            result = *this;
            target = &result;
        }

        int64 pivot_row = 0;
        int64 num_rows = target->rows();
        int64 num_cols = target->cols();

        for (int64 col = 0; col < num_cols; ++col)
        {
            int64 pivot_idx = -1;
            TN max_val = TN(0);

            for (int64 row = pivot_row; row < num_rows; ++row)
            {
                if (std::abs((*target)(row, col)) > std::abs(max_val))
                {
                    max_val = (*target)(row, col);
                    pivot_idx = row;
                }
            }

            if (pivot_idx == -1 || max_val == TN(0)) continue;

            if (pivot_idx != pivot_row)
            {
                target->swap_rows(pivot_row, pivot_idx, true);
            }

            for (int row = pivot_row + 1; row < num_rows; ++row)
            {
                TN factor = (*target)(row, col) / (*target)(pivot_row, col);
                for (int64 j = col; j < num_cols; ++j)
                {
                    (*target)(row, j) -= factor * (*target)(pivot_row, j);
                }
            }

            pivot_row++;
        }
        
        return *target;
    }

    // Constructors

    template <typename TN>
    Matrix<TN>::Matrix() : dense_(), sparse_(), QH_(), betas_() {}

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
        if (i < 0) i += rows();
        if (j < 0) j += cols();
        return dense_(i, j);
    }

    template <typename TN>
    const TN& Matrix<TN>::operator()(int64 i, int64 j) const {
        if (i < 0) i += rows();
        if (j < 0) j += cols();
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
        DenseMatrix<TN> result_dense = dense_ * other.dense_;
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
    Matrix<TN> Matrix<TN>::submatrix(const std::initializer_list<int64>& rows_range,
                                     const std::initializer_list<int64>& cols_range) {
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