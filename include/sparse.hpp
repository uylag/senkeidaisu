#pragma once
#include <dense.hpp>
#include <declarations.hpp>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <algorithm>

namespace senkeidaisu {
template <typename TN>
class DenseMatrix;

template <typename TN>
class SparseMatrix 
{
public:
    friend class DenseMatrix<TN>;

    SparseMatrix();

    SparseMatrix(const int64& rows, const int64& cols, 
                 const TN& default_value = TN(0));

    SparseMatrix(const TN*& c_matrix, const int64& rows, 
                 const int64& cols);

    SparseMatrix(const std::vector<TN>& c_matrix, const int64& rows, 
                 const int64& cols);

    SparseMatrix(const std::vector<std::vector<TN>>& matrix);

    SparseMatrix(const SparseMatrix<TN>& other);

    SparseMatrix(SparseMatrix<TN>&& other) noexcept;

    SparseMatrix(std::initializer_list<std::initializer_list<TN>> matrix);

    SparseMatrix(std::initializer_list<TN> matrix, int rows, int columns);

    SparseMatrix(const DenseMatrix<TN>& other);

    ~SparseMatrix();

    TN& operator()(int64 i, int64 j);
    const TN& operator()(int64 i, int64 j) const;

    SparseMatrix<TN> operator+(const TN& scalar) const;
    SparseMatrix<TN> operator+(const TN*& c_matrix) const;
    SparseMatrix<TN> operator+(const TN** matrix) const;
    SparseMatrix<TN> operator+(const std::vector<TN> c_matrix) const;
    SparseMatrix<TN> operator+(const std::vector<std::vector<TN>> c_matrix) const;
    SparseMatrix<TN> operator+(const SparseMatrix<TN>& other) const;

    SparseMatrix<TN> operator-(const TN& scalar) const;
    SparseMatrix<TN> operator-(const TN*& c_matrix) const;
    SparseMatrix<TN> operator-(const TN** matrix) const;
    SparseMatrix<TN> operator-(const std::vector<TN> c_matrix) const;
    SparseMatrix<TN> operator-(const std::vector<std::vector<TN>> c_matrix) const;
    SparseMatrix<TN> operator-(const SparseMatrix<TN>& other) const;

    SparseMatrix<TN> operator*(const TN& scalar) const;
    SparseMatrix<TN> operator*(const TN*& c_matrix) const;
    SparseMatrix<TN> operator*(const TN** matrix) const;
    SparseMatrix<TN> operator*(const std::vector<TN> c_matrix) const;
    SparseMatrix<TN> operator*(const std::vector<std::vector<TN>> c_matrix) const;
    SparseMatrix<TN> operator*(const SparseMatrix<TN>& other) const;

    SparseMatrix<TN>& operator=(const SparseMatrix<TN>& other);
    SparseMatrix<TN>& operator=(SparseMatrix<TN>&& other) noexcept;

    SparseMatrix<TN> operator-() const;

    SparseMatrix<TN>& operator+=(const TN& scalar);
    SparseMatrix<TN>& operator+=(const TN*& c_matrix);
    SparseMatrix<TN>& operator+=(const TN** matrix);
    SparseMatrix<TN>& operator+=(const std::vector<TN> c_matrix);
    SparseMatrix<TN>& operator+=(const std::vector<std::vector<TN>> c_matrix);
    SparseMatrix<TN>& operator+=(const SparseMatrix<TN>& other);

    SparseMatrix<TN>& operator-=(const TN& scalar);
    SparseMatrix<TN>& operator-=(const TN*& c_matrix);
    SparseMatrix<TN>& operator-=(const TN** matrix);
    SparseMatrix<TN>& operator-=(const std::vector<TN> c_matrix);
    SparseMatrix<TN>& operator-=(const std::vector<std::vector<TN>> c_matrix);
    SparseMatrix<TN>& operator-=(const SparseMatrix<TN>& other);

    SparseMatrix<TN>& operator*=(const TN& scalar);
    SparseMatrix<TN>& operator*=(const TN*& c_matrix);
    SparseMatrix<TN>& operator*=(const TN** matrix);
    SparseMatrix<TN>& operator*=(const std::vector<TN> c_matrix);
    SparseMatrix<TN>& operator*=(const std::vector<std::vector<TN>> c_matrix);
    SparseMatrix<TN>& operator*=(const SparseMatrix<TN>& other);

    bool operator==(const SparseMatrix<TN>& other) const;
    bool operator!=(const SparseMatrix<TN>& other) const { return !(*this == other); };

    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix<TN>& matrix)
    {
        for (int64 i = 0; i < matrix.rows_; ++i)
        {
            for (int64 j = 0; j < matrix.columns_; ++j)
            {
                os << matrix(i, j) << " ";
            }
            os << '\n';
        }
        return os;
    };

    DenseMatrix<TN>                     to_dense() const;

    // Function for turning c_matrix to 2d
    inline std::vector<std::vector<TN>> dd() const;
    inline std::vector<TN>              v() const;

    inline SparseMatrix<TN>              T(bool inplace = false);

    inline bool is_empty() const 
        { return rows_ * columns_ == 0; };

    inline bool is_square() const 
        { return rows_ == columns_; };

    inline int64                        rows() const { return rows_; };
    inline static int64                 rows(const SparseMatrix& matrix) { return matrix.rows_; };

    inline int64                        cols() const { return columns_; };
    inline static int64                 cols(const SparseMatrix& matrix) { return matrix.columns_; };

    inline int64                        size() const { return rows_ * columns_; };
    inline static int64                 size(const SparseMatrix& matrix) { return matrix.size_; };

    inline static SparseMatrix<TN>       identity(int64 size);
    inline static SparseMatrix<TN>       diagonal(const std::vector<TN>& values);

    inline SparseMatrix<TN>              submatrix(const std::initializer_list<int>& rows_start_end,
                                                   const std::initializer_list<int>& columns_start_end);

    inline void                         fill(const TN& value, bool remain_filled = true);
    inline void                         resize(int64 rows, int64 cols);

    SparseMatrix<TN>                    clone() const { return *this; };
    

private:
    static void freeArrays(SparseMatrix<TN>* mat) {
        if(mat->crs_row_ptr_) { delete[] mat->crs_row_ptr_; mat->crs_row_ptr_ = nullptr; }
        if(mat->crs_col_idx_) { delete[] mat->crs_col_idx_; mat->crs_col_idx_ = nullptr; }
        if(mat->crs_values_)  { delete[] mat->crs_values_;  mat->crs_values_  = nullptr; }
        if(mat->ccs_col_ptr_) { delete[] mat->ccs_col_ptr_; mat->ccs_col_ptr_ = nullptr; }
        if(mat->ccs_row_idx_) { delete[] mat->ccs_row_idx_; mat->ccs_row_idx_ = nullptr; }
        if(mat->ccs_values_)  { delete[] mat->ccs_values_;  mat->ccs_values_  = nullptr; }
    }

    int64 rows_, columns_;
    mutable int64* crs_row_ptr_;
    mutable int64* crs_col_idx_;
    mutable TN* crs_values_;

    mutable int64* ccs_col_ptr_;
    mutable int64* ccs_row_idx_;
    mutable TN* ccs_values_;
};

template <typename TN>
DenseMatrix<TN> SparseMatrix<TN>::to_dense() const 
{
    int64 nrows = rows_;
    int64 nnz = crs_row_ptr_[rows_];
    int64 ncols = columns_;
    std::vector<std::vector<TN>> dense(nrows, std::vector<TN>(ncols, 0));
    for (int64 i = 0; i < nrows; ++i) {
        for (int64 idx = crs_row_ptr_[i]; idx < crs_row_ptr_[i+1]; ++idx) {
            int64 col = crs_col_idx_[idx];
            dense[i][col] = crs_values_[idx];
        }
    }
    return DenseMatrix<TN>(dense);
}

template <typename TN>
SparseMatrix<TN>::SparseMatrix()
    : rows_(0), columns_(0),
      crs_row_ptr_(nullptr), crs_col_idx_(nullptr), crs_values_(nullptr),
      ccs_col_ptr_(nullptr), ccs_row_idx_(nullptr), ccs_values_(nullptr)
      {}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(const int64& rows, const int64& cols, const TN& default_value)
    : rows_(rows), columns_(cols)
{
    if (rows < 0 || cols < 0)
        throw std::invalid_argument("Rows and columns must be non-negative");

    int64 total = rows * cols;
    bool nonzero = (default_value != 0);
    int64 nnz = nonzero ? total : 0;

    crs_row_ptr_ = new int64[rows_ + 1];
    for (int64 i = 0; i < rows_; i++) {
        crs_row_ptr_[i] = nonzero ? i * cols : 0;
    }
    crs_row_ptr_[rows_] = nonzero ? total : 0;

    if (nonzero) {
        crs_col_idx_ = new int64[total];
        crs_values_  = new TN[total];
        for (int64 i = 0; i < rows_; i++) {
            for (int64 j = 0; j < cols; j++) {
                crs_col_idx_[i * cols + j] = j;
                crs_values_[i * cols + j]  = default_value;
            }
        }
    } else {
        crs_col_idx_ = nullptr;
        crs_values_  = nullptr;
    }

    ccs_col_ptr_ = new int64[columns_ + 1];
    if (nonzero) {
        ccs_row_idx_ = new int64[total];
        ccs_values_  = new TN[total];
        for (int64 j = 0; j < columns_; j++) {
            ccs_col_ptr_[j] = j * rows_;
            for (int64 i = 0; i < rows_; i++) {
                ccs_row_idx_[j * rows_ + i] = i;
                ccs_values_[j * rows_ + i]  = default_value;
            }
        }
        ccs_col_ptr_[columns_] = total;
    } else {
        ccs_row_idx_ = nullptr;
        ccs_values_  = nullptr;
        for (int64 j = 0; j < columns_ + 1; j++) {
            ccs_col_ptr_[j] = 0;
        }
    }
}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(const TN*& c_matrix, const int64& rows, const int64& cols)
    : rows_(rows), columns_(cols)
{
    if (rows < 0 || cols < 0)
        throw std::invalid_argument("Rows and columns must be non-negative");
    int64 total = rows * cols;
    int64 count = 0;
    for (int64 i = 0; i < total; i++) {
        if (c_matrix[i] != 0)
            count++;
    }
    crs_row_ptr_ = new int64[rows_ + 1];
    crs_col_idx_ = new int64[count];
    crs_values_  = new TN[count];
    int64 pos = 0;
    for (int64 i = 0; i < rows_; i++) {
        crs_row_ptr_[i] = pos;
        for (int64 j = 0; j < cols; j++) {
            int64 index = i * cols + j;
            if (c_matrix[index] != 0) {
                crs_col_idx_[pos] = j;
                crs_values_[pos]  = c_matrix[index];
                pos++;
            }
        }
    }
    crs_row_ptr_[rows_] = pos;
    ccs_col_ptr_ = new int64[columns_ + 1];
    ccs_row_idx_ = new int64[count];
    ccs_values_  = new TN[count];
    std::vector<int64> colCounts(cols, 0);
    for (int64 i = 0; i < count; i++) {
        int64 col = static_cast<int64>(crs_col_idx_[i]);
        colCounts[col]++;
    }
    ccs_col_ptr_[0] = 0;
    for (int64 j = 0; j < cols; j++) {
        ccs_col_ptr_[j+1] = ccs_col_ptr_[j] + colCounts[j];
        colCounts[j] = static_cast<int64>(ccs_col_ptr_[j]);
    }
    for (int64 i = 0; i < rows_; i++) {
        for (int64 k = static_cast<int64>(crs_row_ptr_[i]); k < static_cast<int64>(crs_row_ptr_[i+1]); k++) {
            int64 j = static_cast<int64>(crs_col_idx_[k]);
            int64 pos_ccs = colCounts[j]++;
            ccs_row_idx_[pos_ccs] = i;
            ccs_values_[pos_ccs]  = crs_values_[k];
        }
    }
}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(const std::vector<TN>& c_matrix, const int64& rows, const int64& cols)
    : SparseMatrix(c_matrix.data(), rows, cols) {}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(const std::vector<std::vector<TN>>& matrix)
    : rows_(matrix.size()), columns_(matrix.empty() ? 0 : matrix[0].size())
{
    int64 count = 0;
    for (const auto& row : matrix) {
        if (row.size() != static_cast<size_t>(columns_))
            throw std::invalid_argument("All rows must have the same number of columns");
        for (auto val : row)
            if (val != 0)
                count++;
    }
    crs_row_ptr_ = new int64[rows_ + 1];
    crs_col_idx_ = new int64[count];
    crs_values_  = new TN[count];
    int64 pos = 0;
    for (int64 i = 0; i < rows_; i++) {
        crs_row_ptr_[i] = pos;
        for (int64 j = 0; j < columns_; j++) {
            if (matrix[i][j] != 0) {
                crs_col_idx_[pos] = j;
                crs_values_[pos]  = matrix[i][j];
                pos++;
            }
        }
    }
    crs_row_ptr_[rows_] = pos;

    ccs_col_ptr_ = new int64[columns_ + 1];
    ccs_row_idx_ = new int64[count];
    ccs_values_  = new TN[count];
    std::vector<int64> colCounts(columns_, 0);
    for (int64 i = 0; i < pos; i++) {
        int64 col = static_cast<int64>(crs_col_idx_[i]);
        colCounts[col]++;
    }
    ccs_col_ptr_[0] = 0;
    for (int64 j = 0; j < columns_; j++) {
        ccs_col_ptr_[j+1] = ccs_col_ptr_[j] + colCounts[j];
        colCounts[j] = static_cast<int64>(ccs_col_ptr_[j]);
    }
    for (int64 i = 0; i < rows_; i++) {
        for (int64 k = static_cast<int64>(crs_row_ptr_[i]); k < static_cast<int64>(crs_row_ptr_[i+1]); k++) {
            int64 j = static_cast<int64>(crs_col_idx_[k]);
            int64 pos_ccs = colCounts[j]++;
            ccs_row_idx_[pos_ccs] = i;
            ccs_values_[pos_ccs]  = crs_values_[k];
        }
    }
}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(const SparseMatrix<TN>& other)
    : rows_(other.rows_), columns_(other.columns_)
{
    int64 nnz = static_cast<int64>(other.crs_row_ptr_[other.rows_]);
    crs_row_ptr_ = new int64[rows_ + 1];
    std::copy(other.crs_row_ptr_, other.crs_row_ptr_ + rows_ + 1, crs_row_ptr_);
    if (nnz > 0) {
        crs_col_idx_ = new int64[nnz];
        crs_values_  = new TN[nnz];
        std::copy(other.crs_col_idx_, other.crs_col_idx_ + nnz, crs_col_idx_);
        std::copy(other.crs_values_, other.crs_values_ + nnz, crs_values_);
    } else {
        crs_col_idx_ = nullptr;
        crs_values_  = nullptr;
    }
    ccs_col_ptr_ = new int64[columns_ + 1];
    std::copy(other.ccs_col_ptr_, other.ccs_col_ptr_ + columns_ + 1, ccs_col_ptr_);
    if (nnz > 0) {
        ccs_row_idx_ = new int64[nnz];
        ccs_values_  = new TN[nnz];
        std::copy(other.ccs_row_idx_, other.ccs_row_idx_ + nnz, ccs_row_idx_);
        std::copy(other.ccs_values_, other.ccs_values_ + nnz, ccs_values_);
    } else {
        ccs_row_idx_ = nullptr;
        ccs_values_  = nullptr;
    }
}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(SparseMatrix<TN>&& other) noexcept
    : rows_(other.rows_), columns_(other.columns_),
      crs_row_ptr_(other.crs_row_ptr_), crs_col_idx_(other.crs_col_idx_), crs_values_(other.crs_values_),
      ccs_col_ptr_(other.ccs_col_ptr_), ccs_row_idx_(other.ccs_row_idx_), ccs_values_(other.ccs_values_)
{
    other.rows_ = 0;
    other.columns_ = 0;
    other.crs_row_ptr_ = nullptr;
    other.crs_col_idx_ = nullptr;
    other.crs_values_  = nullptr;
    other.ccs_col_ptr_ = nullptr;
    other.ccs_row_idx_ = nullptr;
    other.ccs_values_  = nullptr;
}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(std::initializer_list<std::initializer_list<TN>> matrix)
{
    rows_ = matrix.size();
    columns_ = matrix.begin()->size();
    int64 count = 0;
    for (auto& row : matrix) {
        if (row.size() != columns_)
            throw std::invalid_argument("All rows must have the same number of columns");
        for (auto val : row)
            if (val != 0)
                count++;
    }
    crs_row_ptr_ = new int64[rows_ + 1];
    crs_col_idx_ = new int64[count];
    crs_values_  = new TN[count];
    int64 pos = 0;
    int64 i = 0;
    for (auto& row : matrix) {
        crs_row_ptr_[i] = pos;
        int64 j = 0;
        for (auto val : row) {
            if (val != 0) {
                crs_col_idx_[pos] = j;
                crs_values_[pos]  = val;
                pos++;
            }
           	j++;
        }
        i++;
    }
    crs_row_ptr_[rows_] = pos;

    ccs_col_ptr_ = new int64[columns_ + 1];
    ccs_row_idx_ = new int64[count];
    ccs_values_  = new TN[count];
    std::vector<int64> colCounts(columns_, 0);
    for (int64 k = 0; k < pos; k++) {
        int64 col = static_cast<int64>(crs_col_idx_[k]);
        colCounts[col]++;
    }
    ccs_col_ptr_[0] = 0;
    for (int64 j = 0; j < columns_; j++) {
        ccs_col_ptr_[j+1] = ccs_col_ptr_[j] + colCounts[j];
        colCounts[j] = static_cast<int64>(ccs_col_ptr_[j]);
    }
    for (int64 i = 0; i < rows_; i++) {
        for (int64 k = static_cast<int64>(crs_row_ptr_[i]); k < static_cast<int64>(crs_row_ptr_[i+1]); k++) {
            int64 j = static_cast<int64>(crs_col_idx_[k]);
            int64 pos_ccs = colCounts[j]++;
            ccs_row_idx_[pos_ccs] = i;
            ccs_values_[pos_ccs]  = crs_values_[k];
        }
    }
}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(std::initializer_list<TN> matrix, int rows, int columns)
    : rows_(rows), columns_(columns)
{
    if (rows < 0 || columns < 0)
        throw std::invalid_argument("Rows and columns must be non-negative");
    if (static_cast<int64>(matrix.size()) != rows * columns)
        throw std::invalid_argument("Initializer list size does not match matrix dimensions");
    std::vector<TN> vec(matrix);
    *this = SparseMatrix<TN>(vec, rows, columns);
}

template <typename TN>
SparseMatrix<TN>::SparseMatrix(const DenseMatrix<TN>& other)
    : rows_(other.rows()), columns_(other.cols())
{
    int64 total = rows_ * columns_;
    int64 count = 0;
    for (int64 i = 0; i < rows_; i++) {
        for (int64 j = 0; j < columns_; j++) {
            if (other(i, j) != 0)
                count++;
        }
    }
    crs_row_ptr_ = new int64[rows_ + 1];
    crs_col_idx_ = new int64[count];
    crs_values_  = new TN[count];
    int64 pos = 0;
    for (int64 i = 0; i < rows_; i++) {
        crs_row_ptr_[i] = pos;
        for (int64 j = 0; j < columns_; j++) {
            TN val = other(i, j);
            if (val != 0) {
                crs_col_idx_[pos] = j;
                crs_values_[pos]  = val;
                pos++;
            }
        }
    }
    crs_row_ptr_[rows_] = pos;
    ccs_col_ptr_ = new int64[columns_ + 1];
    ccs_row_idx_ = new int64[count];
    ccs_values_  = new TN[count];
    std::vector<int64> colCounts(columns_, 0);
    for (int64 k = 0; k < pos; k++) {
        int64 col = static_cast<int64>(crs_col_idx_[k]);
        colCounts[col]++;
    }
    ccs_col_ptr_[0] = 0;
    for (int64 j = 0; j < columns_; j++) {
        ccs_col_ptr_[j+1] = ccs_col_ptr_[j] + colCounts[j];
        colCounts[j] = static_cast<int64>(ccs_col_ptr_[j]);
    }
    for (int64 i = 0; i < rows_; i++) {
        for (int64 k = static_cast<int64>(crs_row_ptr_[i]); k < static_cast<int64>(crs_row_ptr_[i+1]); k++) {
            int64 j = static_cast<int64>(crs_col_idx_[k]);
            int64 pos_ccs = colCounts[j]++;
            ccs_row_idx_[pos_ccs] = i;
            ccs_values_[pos_ccs]  = crs_values_[k];
        }
    }
}

template <typename TN>
SparseMatrix<TN>::~SparseMatrix() {
    freeArrays(this);
}

template <typename TN>
TN& SparseMatrix<TN>::operator()(int64 i, int64 j) {
    if(i < 0 || i >= rows_ || j < 0 || j >= columns_) throw std::out_of_range("");
    int64 rs = static_cast<int64>(crs_row_ptr_[i]);
    int64 re = static_cast<int64>(crs_row_ptr_[i+1]);
    for (int64 k = rs; k < re; k++) {
        if(static_cast<int64>(crs_col_idx_[k]) == j) return crs_values_[k];
    }
    auto insertElement = [this](int64 i, int64 j, TN val) -> TN& {
        int64 oldNnz = static_cast<int64>(crs_row_ptr_[rows_]);
        int64 pos = static_cast<int64>(crs_row_ptr_[i]);
        while(pos < static_cast<int64>(crs_row_ptr_[i+1]) && static_cast<int64>(crs_col_idx_[pos]) < j) pos++;
        int64* new_crs_row_ptr = new int64[rows_+1];
        int64* new_crs_col_idx = new int64[oldNnz+1];
        TN* new_crs_values = new TN[oldNnz+1];
        new_crs_row_ptr[0] = crs_row_ptr_[0];
        for(int64 r = 0; r < rows_; r++){
            new_crs_row_ptr[r+1] = crs_row_ptr_[r+1] + (r >= i ? 1 : 0);
        }
        for (int64 k = 0; k < pos; k++){
            new_crs_col_idx[k] = crs_col_idx_[k];
            new_crs_values[k] = crs_values_[k];
        }
        new_crs_col_idx[pos] = j;
        new_crs_values[pos] = val;
        for (int64 k = pos; k < oldNnz; k++){
            new_crs_col_idx[k+1] = crs_col_idx_[k];
            new_crs_values[k+1] = crs_values_[k];
        }
        delete[] crs_row_ptr_;
        delete[] crs_col_idx_;
        delete[] crs_values_;
        crs_row_ptr_ = new_crs_row_ptr;
        crs_col_idx_ = new_crs_col_idx;
        crs_values_ = new_crs_values;
        int64 newNnz = oldNnz + 1;
        int64* new_ccs_col_ptr = new int64[columns_+1];
        int64* new_ccs_row_idx = new int64[newNnz];
        TN* new_ccs_values = new TN[newNnz];
        std::vector<int64> cc(columns_,0);
        for(int64 k = 0; k < newNnz; k++){
            int64 col = static_cast<int64>(crs_col_idx_[k]);
            cc[col]++;
        }
        new_ccs_col_ptr[0] = 0;
        for(int64 jcol = 0; jcol < columns_; jcol++){
            new_ccs_col_ptr[jcol+1] = new_ccs_col_ptr[jcol] + cc[jcol];
            cc[jcol] = new_ccs_col_ptr[jcol];
        }
        for(int64 r = 0; r < rows_; r++){
            for(int64 k = static_cast<int64>(crs_row_ptr_[r]); k < static_cast<int64>(crs_row_ptr_[r+1]); k++){
                int64 col = static_cast<int64>(crs_col_idx_[k]);
                int64 pos_ccs = cc[col]++;
                new_ccs_row_idx[pos_ccs] = r;
                new_ccs_values[pos_ccs] = crs_values_[k];
            }
        }
        if(ccs_col_ptr_) { delete[] ccs_col_ptr_; }
        if(ccs_row_idx_) { delete[] ccs_row_idx_; }
        if(ccs_values_) { delete[] ccs_values_; }
        ccs_col_ptr_ = new_ccs_col_ptr;
        ccs_row_idx_ = new_ccs_row_idx;
        ccs_values_ = new_ccs_values;
        return crs_values_[pos];
    };
    return insertElement(i, j, 0);
}

template <typename TN>
const TN& SparseMatrix<TN>::operator()(int64 i, int64 j) const {
    if(i < 0 || i >= rows_ || j < 0 || j >= columns_) throw std::out_of_range("");
    int64 rs = static_cast<int64>(crs_row_ptr_[i]);
    int64 re = static_cast<int64>(crs_row_ptr_[i+1]);
    for (int64 k = rs; k < re; k++) {
        if(static_cast<int64>(crs_col_idx_[k]) == j) return crs_values_[k];
    }
    static TN zero = 0;
    return zero;
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator+(const TN& scalar) const {
    if(scalar == 0) return *this;
    DenseMatrix<TN> d = this->to_dense();
    d = d + scalar;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator+(const TN*& c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix, rows_, columns_);
    d = d + o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator+(const TN** matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(matrix, rows_, columns_);
    d = d + o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator+(const std::vector<TN> c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix, rows_, columns_);
    d = d + o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator+(const std::vector<std::vector<TN>> c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix);
    d = d + o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator+(const SparseMatrix<TN>& other) const {
    if(rows_ != other.rows_ || columns_ != other.columns_) throw std::invalid_argument("");
    int64 rs, re, pos1, pos2;
    std::vector<TN> nc;
    std::vector<TN> nv;
    std::vector<TN> np(rows_+1,0);
    for(int64 i=0;i<rows_;i++){
        pos1 = static_cast<int64>(crs_row_ptr_[i]);
        pos2 = static_cast<int64>(other.crs_row_ptr_[i]);
        re = static_cast<int64>(crs_row_ptr_[i+1]);
        int64 re2 = static_cast<int64>(other.crs_row_ptr_[i+1]);
        while(pos1 < re || pos2 < re2){
            if(pos1 < re && (pos2 >= re2 || static_cast<int64>(crs_col_idx_[pos1]) < static_cast<int64>(other.crs_col_idx_[pos2]))){
                nc.push_back(crs_col_idx_[pos1]);
                nv.push_back(crs_values_[pos1]);
                pos1++;
            } else if(pos2 < re2 && (pos1 >= re || static_cast<int64>(other.crs_col_idx_[pos2]) < static_cast<int64>(crs_col_idx_[pos1]))){
                nc.push_back(other.crs_col_idx_[pos2]);
                nv.push_back(other.crs_values_[pos2]);
                pos2++;
            } else {
                TN sum = crs_values_[pos1] + other.crs_values_[pos2];
                if(sum != 0){
                    nc.push_back(crs_col_idx_[pos1]);
                    nv.push_back(sum);
                }
                pos1++; pos2++;
            }
        }
        np[i+1] = nc.size();
    }
    SparseMatrix<TN> r;
    r.rows_ = rows_;
    r.columns_ = columns_;
    r.crs_row_ptr_ = new int64[rows_+1];
    r.crs_col_idx_ = new int64[nc.size()];
    r.crs_values_ = new TN[nv.size()];
    std::copy(np.begin(), np.end(), r.crs_row_ptr_);
    std::copy(nc.begin(), nc.end(), r.crs_col_idx_);
    std::copy(nv.begin(), nv.end(), r.crs_values_);
    int64 total = r.crs_row_ptr_[rows_];
    r.ccs_col_ptr_ = new int64[columns_+1];
    r.ccs_row_idx_ = new int64[total];
    r.ccs_values_ = new TN[total];
    std::vector<int64> cc(columns_,0);
    for(int64 k=0;k<total;k++){
        int64 col = static_cast<int64>(r.crs_col_idx_[k]);
        cc[col]++;
    }
    r.ccs_col_ptr_[0] = 0;
    for(int64 j=0;j<columns_;j++){
        r.ccs_col_ptr_[j+1] = r.ccs_col_ptr_[j] + cc[j];
        cc[j] = r.ccs_col_ptr_[j];
    }
    for(int64 i=0;i<rows_;i++){
        for(int64 k = static_cast<int64>(r.crs_row_ptr_[i]); k < static_cast<int64>(r.crs_row_ptr_[i+1]); k++){
            int64 col = static_cast<int64>(r.crs_col_idx_[k]);
            int64 pos_ccs = cc[col]++;
            r.ccs_row_idx_[pos_ccs] = i;
            r.ccs_values_[pos_ccs] = r.crs_values_[k];
        }
    }
    return r;
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator-(const TN& scalar) const {
    return (*this) + (-scalar);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator-(const TN*& c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix, rows_, columns_);
    d = d - o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator-(const TN** matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(matrix, rows_, columns_);
    d = d - o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator-(const std::vector<TN> c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix, rows_, columns_);
    d = d - o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator-(const std::vector<std::vector<TN>> c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix);
    d = d - o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator-(const SparseMatrix<TN>& other) const {
    return (*this) + (-other);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator*(const TN& scalar) const {
    SparseMatrix<TN> r(*this);
    if(scalar == 0){
        int64 total = static_cast<int64>(r.crs_row_ptr_[r.rows_]);
        delete[] r.crs_col_idx_;
        delete[] r.crs_values_;
        r.crs_col_idx_ = nullptr;
        r.crs_values_ = nullptr;
        for(int64 i=0;i<=r.rows_;i++){
            r.crs_row_ptr_[i] = 0;
        }
        delete[] r.ccs_col_ptr_;
        delete[] r.ccs_row_idx_;
        delete[] r.ccs_values_;
        r.ccs_col_ptr_ = new int64[r.columns_+1];
        r.ccs_row_idx_ = nullptr;
        r.ccs_values_ = nullptr;
        for(int64 j=0;j<=r.columns_;j++){
            r.ccs_col_ptr_[j] = 0;
        }
    } else {
        int64 total = static_cast<int64>(r.crs_row_ptr_[r.rows_]);
        for(int64 k=0;k<total;k++){
            r.crs_values_[k] *= scalar;
        }
        total = static_cast<int64>(r.ccs_col_ptr_[r.columns_]);
        for(int64 k=0;k<total;k++){
            r.ccs_values_[k] *= scalar;
        }
    }
    return r;
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator*(const TN*& c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix, rows_, columns_);
    d = d * o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator*(const TN** matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(matrix, rows_, columns_);
    d = d * o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator*(const std::vector<TN> c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix, rows_, columns_);
    d = d * o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator*(const std::vector<std::vector<TN>> c_matrix) const {
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> o(c_matrix);
    d = d * o;
    return SparseMatrix<TN>(d);
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator*(const SparseMatrix<TN>& other) const {
    if (columns_ != other.rows_) throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    
    SparseMatrix<TN> result(rows_, other.columns_);
    std::vector<std::vector<std::pair<int64, TN>>> temp_result(rows_); 

    for (int64 i = 0; i < rows_; ++i) {
        for (int64 k = crs_row_ptr_[i]; k < crs_row_ptr_[i + 1]; ++k) {
            int64 col = crs_col_idx_[k]; 
            TN val = crs_values_[k];     

            for (int64 j = 0; j < other.columns_; ++j) {
                for (int64 p = other.ccs_col_ptr_[j]; p < other.ccs_col_ptr_[j + 1]; ++p) {
                    if (other.ccs_row_idx_[p] == col) { 
                        TN product = val * other.ccs_values_[p];
                        temp_result[i].push_back({j, product});
                    }
                }
            }
        }
    }

    std::vector<int64> nc; 
    std::vector<TN> nv;    
    std::vector<int64> np(rows_ + 1, 0); // crs_row_ptr

    for (int64 i = 0; i < rows_; ++i) {
        std::sort(temp_result[i].begin(), temp_result[i].end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        int64 last_j = -1;
        TN sum = 0;
        for (const auto& [j, val] : temp_result[i]) {
            if (j != last_j && last_j != -1 && sum != 0) {
                nc.push_back(last_j);
                nv.push_back(sum);
            }
            if (j != last_j) {
                sum = val;
                last_j = j;
            } else {
                sum += val;
            }
        }
        if (last_j != -1 && sum != 0) {
            nc.push_back(last_j);
            nv.push_back(sum);
        }
        np[i + 1] = nc.size();
    }

    result.crs_row_ptr_ = new int64[rows_ + 1];
    result.crs_col_idx_ = new int64[nc.size()];
    result.crs_values_ = new TN[nv.size()];
    std::copy(np.begin(), np.end(), result.crs_row_ptr_);
    std::copy(nc.begin(), nc.end(), result.crs_col_idx_);
    std::copy(nv.begin(), nv.end(), result.crs_values_);

    int64 nnz = nc.size();
    result.ccs_col_ptr_ = new int64[other.columns_ + 1];
    result.ccs_row_idx_ = new int64[nnz];
    result.ccs_values_ = new TN[nnz];
    std::vector<int64> col_counts(other.columns_, 0);
    for (int64 k = 0; k < nnz; ++k) {
        col_counts[result.crs_col_idx_[k]]++;
    }
    result.ccs_col_ptr_[0] = 0;
    for (int64 j = 0; j < other.columns_; ++j) {
        result.ccs_col_ptr_[j + 1] = result.ccs_col_ptr_[j] + col_counts[j];
        col_counts[j] = result.ccs_col_ptr_[j];
    }
    for (int64 i = 0; i < rows_; ++i) {
        for (int64 k = result.crs_row_ptr_[i]; k < result.crs_row_ptr_[i + 1]; ++k) {
            int64 j = result.crs_col_idx_[k];
            int64 pos = col_counts[j]++;
            result.ccs_row_idx_[pos] = i;
            result.ccs_values_[pos] = result.crs_values_[k];
        }
    }

    return result;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator=(const SparseMatrix<TN>& other) {
    if(this != &other){
        rows_ = other.rows_;
        columns_ = other.columns_;
        int64 nnz = static_cast<int64>(other.crs_row_ptr_[other.rows_]);
        if(crs_row_ptr_) delete[] crs_row_ptr_;
        crs_row_ptr_ = new int64[rows_+1];
        std::copy(other.crs_row_ptr_, other.crs_row_ptr_+rows_+1, crs_row_ptr_);
        if(nnz > 0){
            if(crs_col_idx_) delete[] crs_col_idx_;
            if(crs_values_) delete[] crs_values_;
            crs_col_idx_ = new int64[nnz];
            crs_values_ = new TN[nnz];
            std::copy(other.crs_col_idx_, other.crs_col_idx_+nnz, crs_col_idx_);
            std::copy(other.crs_values_, other.crs_values_+nnz, crs_values_);
        } else {
            crs_col_idx_ = nullptr;
            crs_values_ = nullptr;
        }
        if(ccs_col_ptr_) delete[] ccs_col_ptr_;
        ccs_col_ptr_ = new int64[columns_+1];
        std::copy(other.ccs_col_ptr_, other.ccs_col_ptr_+columns_+1, ccs_col_ptr_);
        if(nnz > 0){
            if(ccs_row_idx_) delete[] ccs_row_idx_;
            if(ccs_values_) delete[] ccs_values_;
            ccs_row_idx_ = new int64[nnz];
            ccs_values_ = new TN[nnz];
            std::copy(other.ccs_row_idx_, other.ccs_row_idx_+nnz, ccs_row_idx_);
            std::copy(other.ccs_values_, other.ccs_values_+nnz, ccs_values_);
        } else {
            ccs_row_idx_ = nullptr;
            ccs_values_ = nullptr;
        }
    }
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator=(SparseMatrix<TN>&& other) noexcept {
    if(this != &other){
        rows_ = other.rows_;
        columns_ = other.columns_;
        crs_row_ptr_ = other.crs_row_ptr_;
        crs_col_idx_ = other.crs_col_idx_;
        crs_values_ = other.crs_values_;
        ccs_col_ptr_ = other.ccs_col_ptr_;
        ccs_row_idx_ = other.ccs_row_idx_;
        ccs_values_ = other.ccs_values_;
        other.rows_ = 0;
        other.columns_ = 0;
        other.crs_row_ptr_ = nullptr;
        other.crs_col_idx_ = nullptr;
        other.crs_values_ = nullptr;
        other.ccs_col_ptr_ = nullptr;
        other.ccs_row_idx_ = nullptr;
        other.ccs_values_ = nullptr;
    }
    return *this;
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::operator-() const {
    SparseMatrix<TN> r(*this);
    int64 total = static_cast<int64>(r.crs_row_ptr_[r.rows_]);
    for(int64 k=0;k<total;k++){
        r.crs_values_[k] = -r.crs_values_[k];
    }
    total = static_cast<int64>(r.ccs_col_ptr_[r.columns_]);
    for(int64 k=0;k<total;k++){
        r.ccs_values_[k] = -r.ccs_values_[k];
    }
    return r;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator+=(const TN& scalar) {
    *this = *this + scalar;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator+=(const TN*& c_matrix) {
    *this = *this + c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator+=(const TN** matrix) {
    *this = *this + matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator+=(const std::vector<TN> c_matrix) {
    *this = *this + c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator+=(const std::vector<std::vector<TN>> c_matrix) {
    *this = *this + c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator+=(const SparseMatrix<TN>& other) {
    *this = *this + other;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator-=(const TN& scalar) {
    *this = *this - scalar;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator-=(const TN*& c_matrix) {
    *this = *this - c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator-=(const TN** matrix) {
    *this = *this - matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator-=(const std::vector<TN> c_matrix) {
    *this = *this - c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator-=(const std::vector<std::vector<TN>> c_matrix) {
    *this = *this - c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator-=(const SparseMatrix<TN>& other) {
    *this = *this - other;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator*=(const TN& scalar) {
    *this = *this * scalar;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator*=(const TN*& c_matrix) {
    *this = *this * c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator*=(const TN** matrix) {
    *this = *this * matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator*=(const std::vector<TN> c_matrix) {
    *this = *this * c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator*=(const std::vector<std::vector<TN>> c_matrix) {
    *this = *this * c_matrix;
    return *this;
}

template <typename TN>
SparseMatrix<TN>& SparseMatrix<TN>::operator*=(const SparseMatrix<TN>& other) {
    *this = *this * other;
    return *this;
}

template <typename TN>
bool SparseMatrix<TN>::operator==(const SparseMatrix<TN>& other) const {
    if(rows_ != other.rows_ || columns_ != other.columns_) return false;
    int64 nnz1 = static_cast<int64>(crs_row_ptr_[rows_]);
    int64 nnz2 = static_cast<int64>(other.crs_row_ptr_[other.rows_]);
    if(nnz1 != nnz2) return false;
    for(int64 k=0;k<nnz1;k++){
        if(crs_col_idx_[k] != other.crs_col_idx_[k] || crs_values_[k] != other.crs_values_[k])
            return false;
    }
    return true;
}

template <typename TN>
std::vector<std::vector<TN>> SparseMatrix<TN>::dd() const {
    DenseMatrix<TN> d = this->to_dense();
    return d.dd();
}

template <typename TN>
std::vector<TN> SparseMatrix<TN>::v() const {
    DenseMatrix<TN> d = this->to_dense();
    return d.v();
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::T(bool inplace) {
    SparseMatrix<TN> t;
    int64 tot = static_cast<int64>(ccs_col_ptr_[columns_]);
    int64* n_crs_row = new int64[columns_+1];
    int64* n_crs_col = new int64[tot];
    TN* n_crs_val = new TN[tot];
    n_crs_row[0] = 0;
    for(int64 j=0;j<columns_;j++){
        int64 cnt = static_cast<int64>(ccs_col_ptr_[j+1] - ccs_col_ptr_[j]);
        n_crs_row[j+1] = n_crs_row[j] + cnt;
    }
    for(int64 j=0;j<columns_;j++){
        for(int64 k = static_cast<int64>(ccs_col_ptr_[j]); k < static_cast<int64>(ccs_col_ptr_[j+1]); k++){
            int64 pos = static_cast<int64>(ccs_col_ptr_[j]) + (k - ccs_col_ptr_[j]);
            n_crs_col[pos] = ccs_row_idx_[k];
            n_crs_val[pos] = ccs_values_[k];
        }
    }
    if(inplace){
        delete[] crs_row_ptr_;
        delete[] crs_col_idx_;
        delete[] crs_values_;
        crs_row_ptr_ = n_crs_row;
        crs_col_idx_ = n_crs_col;
        crs_values_ = n_crs_val;
        int64 tmp = rows_;
        rows_ = columns_;
        columns_ = tmp;
        if(ccs_col_ptr_){ delete[] ccs_col_ptr_; ccs_col_ptr_ = nullptr; }
        if(ccs_row_idx_){ delete[] ccs_row_idx_; ccs_row_idx_ = nullptr; }
        if(ccs_values_){ delete[] ccs_values_; ccs_values_ = nullptr; }
        return *this;
    } else {
        t.rows_ = columns_;
        t.columns_ = rows_;
        t.crs_row_ptr_ = n_crs_row;
        t.crs_col_idx_ = n_crs_col;
        t.crs_values_ = n_crs_val;
        int64 tot2 = n_crs_row[columns_];
        t.ccs_col_ptr_ = new int64[rows_+1];
        t.ccs_row_idx_ = new int64[tot2];
        t.ccs_values_ = new TN[tot2];
        std::vector<int64> cc(rows_,0);
        for(int64 k=0;k<tot2;k++){
            int64 col = static_cast<int64>(t.crs_col_idx_[k]);
            cc[col]++;
        }
        t.ccs_col_ptr_[0] = 0;
        for(int64 j=0;j<rows_;j++){
            t.ccs_col_ptr_[j+1] = t.ccs_col_ptr_[j] + cc[j];
            cc[j] = t.ccs_col_ptr_[j];
        }
        for(int64 i=0;i<t.rows_;i++){
            for(int64 k = static_cast<int64>(t.crs_row_ptr_[i]); k < static_cast<int64>(t.crs_row_ptr_[i+1]); k++){
                int64 col = static_cast<int64>(t.crs_col_idx_[k]);
                int64 pos_ccs = cc[col]++;
                t.ccs_row_idx_[pos_ccs] = i;
                t.ccs_values_[pos_ccs] = t.crs_values_[k];
            }
        }
        return t;
    }
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::identity(int64 size) {
    SparseMatrix<TN> r(size, size, 0);
    for(int64 i=0;i<size;i++){
        r(i,i) = 1;
    }
    return r;
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::diagonal(const std::vector<TN>& values) {
    int64 size = values.size();
    SparseMatrix<TN> r(size, size, 0);
    for(int64 i=0;i<size;i++){
        r(i,i) = values[i];
    }
    return r;
}

template <typename TN>
SparseMatrix<TN> SparseMatrix<TN>::submatrix(const std::initializer_list<int>& rs_re, 
                                             const std::initializer_list<int>& cs_ce) {
    
    if (rs_re.size() != 2 || cs_ce.size() != 2) {
        throw std::invalid_argument("Ranges must have exactly 2 elements: start and end.");
    }

    auto rs = *rs_re.begin();
    auto re = *(rs_re.begin() + 1);
    auto cs = *cs_ce.begin();
    auto ce = *(cs_ce.begin() + 1);
    
    if(rs < 0 || rs > re || re > rows_ || cs < 0 || cs > ce || ce > columns_) 
        throw std::out_of_range("Problem in submatrix creation: invalid rs_re or cs_ce");
    DenseMatrix<TN> d = this->to_dense();
    DenseMatrix<TN> s = d.submatrix({rs,re}, {cs,ce});
    return SparseMatrix<TN>(s);
}

template <typename TN>
void SparseMatrix<TN>::fill(const TN& value, bool rf) {
    DenseMatrix<TN> d = this->to_dense();
    d.fill(value, rf);
    *this = SparseMatrix<TN>(d);
}

template <typename TN>
void SparseMatrix<TN>::resize(int64 r, int64 c) {
    DenseMatrix<TN> d = this->to_dense();
    d.resize(r,c);
    *this = SparseMatrix<TN>(d);
}
}

