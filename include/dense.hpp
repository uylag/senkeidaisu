#pragma once
#include <declarations.hpp>
#include <sparse.hpp>
#include <algorithm>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <utility>

namespace senkeidaisu 
{
    template <typename TN>
    class SparseMatrix;

    template <typename TN>
    class DenseMatrix 
    {
    public:

        // CONSTRUCTORS & DESTRUCTOR BLOCK
        DenseMatrix();

        DenseMatrix(const int64& rows, const int64& columns, 
                    const TN& default_value = TN(0));

        DenseMatrix(const TN* c_matrix, 
                    const int64& rows, const int64& columns);

        DenseMatrix(const std::vector<TN> c_matrix,
                    const int64& rows, const int64& columns);

        DenseMatrix(const std::vector<std::vector<TN>>& matrix);

        DenseMatrix(const DenseMatrix<TN>& other);

        DenseMatrix(DenseMatrix<TN>&& other) noexcept;

        DenseMatrix(std::initializer_list<std::initializer_list<TN>> matrix);

        DenseMatrix(std::initializer_list<TN> matrix, int rows, int columns);

        DenseMatrix(const SparseMatrix<TN>& other);

        ~DenseMatrix();

        // OPERATIONS

        TN& operator()(int64 i, int64 j);
        const TN& operator()(int64 i, int64 j) const;

        DenseMatrix<TN> operator+(const TN& scalar) const;
        DenseMatrix<TN> operator+(const TN*& c_matrix) const;
        DenseMatrix<TN> operator+(const TN** matrix) const;
        DenseMatrix<TN> operator+(const std::vector<TN> c_matrix) const;
        DenseMatrix<TN> operator+(const std::vector<std::vector<TN>> c_matrix) const;
        DenseMatrix<TN> operator+(const DenseMatrix<TN>& other) const;

        DenseMatrix<TN> operator-(const TN& scalar) const;
        DenseMatrix<TN> operator-(const TN*& c_matrix) const;
        DenseMatrix<TN> operator-(const TN** matrix) const;
        DenseMatrix<TN> operator-(const std::vector<TN> c_matrix) const;
        DenseMatrix<TN> operator-(const std::vector<std::vector<TN>> c_matrix) const;
        DenseMatrix<TN> operator-(const DenseMatrix<TN>& other) const;

        DenseMatrix<TN> operator*(const TN& scalar) const;
        DenseMatrix<TN> operator*(const TN*& c_matrix) const;
        DenseMatrix<TN> operator*(const TN** matrix) const;
        DenseMatrix<TN> operator*(const std::vector<TN> c_matrix) const;
        DenseMatrix<TN> operator*(const std::vector<std::vector<TN>> c_matrix) const;
        DenseMatrix<TN> operator*(const DenseMatrix<TN>& other) const;

        DenseMatrix<TN>& operator=(const DenseMatrix<TN>& other);
        DenseMatrix<TN>& operator=(DenseMatrix<TN>&& other) noexcept;

        DenseMatrix<TN> operator-() const;

        DenseMatrix<TN>& operator+=(const TN& scalar);
        DenseMatrix<TN>& operator+=(const TN*& c_matrix);
        DenseMatrix<TN>& operator+=(const TN** matrix);
        DenseMatrix<TN>& operator+=(const std::vector<TN> c_matrix);
        DenseMatrix<TN>& operator+=(const std::vector<std::vector<TN>> c_matrix);
        DenseMatrix<TN>& operator+=(const DenseMatrix<TN>& other);

        DenseMatrix<TN>& operator-=(const TN& scalar);
        DenseMatrix<TN>& operator-=(const TN*& c_matrix);
        DenseMatrix<TN>& operator-=(const TN** matrix);
        DenseMatrix<TN>& operator-=(const std::vector<TN> c_matrix);
        DenseMatrix<TN>& operator-=(const std::vector<std::vector<TN>> c_matrix);
        DenseMatrix<TN>& operator-=(const DenseMatrix<TN>& other);

        DenseMatrix<TN>& operator*=(const TN& scalar);
        DenseMatrix<TN>& operator*=(const TN*& c_matrix);
        DenseMatrix<TN>& operator*=(const TN** matrix);
        DenseMatrix<TN>& operator*=(const std::vector<TN> c_matrix);
        DenseMatrix<TN>& operator*=(const std::vector<std::vector<TN>> c_matrix);
        DenseMatrix<TN>& operator*=(const DenseMatrix<TN>& other);

        DenseMatrix<TN> operator/(const TN& scalar) const
        {
            if (scalar == TN(0))
                throw std::invalid_argument("Division by zero\n");

            DenseMatrix<TN> result(*this);
            for (int64 i = 0; i < rows(); ++i)
            {
                for (int64 j = 0; j < cols(); ++j)
                {
                    result.c_matrix_[i * row_stride_ + j] /= scalar;
                }
            }
            return result;
        };
        DenseMatrix<TN>& operator/=(const TN& scalar)
        {
            if (scalar == TN(0)) {
                throw std::invalid_argument("Division by zero");
            }
        
            // for (auto& elem : c_matrix_)  // Directly modify `this`
            // {
            //     elem /= scalar;
            // }

            for (int64 i = 0; i < this->rows(); ++i)
            {
                for (int64 j = 0; j < this->cols(); ++j)
                {
                    c_matrix_[i * row_stride_ + j] /= scalar;
                }
            }

            return *this;
        };

        bool operator==(const DenseMatrix<TN>& other) const;
        bool operator!=(const DenseMatrix<TN>& other) const { return !(*this == other); };

        friend std::ostream& operator<<(std::ostream& os, const DenseMatrix<TN>& matrix)
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

        // FUNCTIONS
        
        // Function for turning c_matrix to 2d
        inline std::vector<std::vector<TN>> dd() const;
        inline std::vector<TN>              v() const;

        inline DenseMatrix<TN>              T(bool inplace = false);

        inline bool is_empty() const 
            { return size_ == 0; };

        inline bool is_square() const 
            { return rows_ == columns_; };

        inline int64                        rows() const { return rows_; };
        inline static int64                 rows(const DenseMatrix& matrix) { return matrix.rows_; };

        inline int64                        cols() const { return columns_; };
        inline static int64                 cols(const DenseMatrix& matrix) { return matrix.cols_; };

        inline int64                        size() const { return size_; };
        inline static int64                 size(const DenseMatrix& matrix) { return matrix.size_; };

        // Scalar multiplication
        inline static DenseMatrix<TN>       identity(int64 size);
        inline static DenseMatrix<TN>       diagonal(const std::vector<TN>& values);

        inline DenseMatrix<TN>              submatrix(const std::pair<int64, int64>& rows_start_end,
                                                      const std::pair<int64, int64>& cols_start_end);

        inline void                         fill(const TN& value);
        inline void                         resize(int64 rows, int64 cols, const TN& def_value = TN(0));

        inline TN*                          begin() { return c_matrix_; };
        inline TN*                          end() { return c_matrix_ + size_; };
        inline const TN*                    begin() const { return c_matrix_; };
        inline const TN*                    end() const { return c_matrix_ + size_; };

        inline TN*                          data() { return c_matrix_; };
        inline const TN*                    data() const { return c_matrix_; };

        DenseMatrix<TN>                     clone() const { return *this; };

    protected:
        
    private:
        // Private functions
        void                                _allocate(int64 capacity);
        void                                _allocate(int64 rows, int64 cols);
        
        // === MATRIX ===
        TN* c_matrix_;

        // === INT64 aka long long instances ===
        int64 rows_, columns_, capacity_, size_, row_stride_;

        // === BOOL INSTANCES ===
        bool memory_owner_;
    };

    // ============================================
    // |                                          |
    // |                                          |
    // |             IMPLEMENTATION               |
    // |                                          |
    // |                                          |
    // ============================================

    template <typename TN>
    void DenseMatrix<TN>::fill(const TN& val)
    {
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                // Use the correct linear index into the one-dimensional array.
                if (c_matrix_[i * row_stride_ + j] == TN(0)) {
                    c_matrix_[i * row_stride_ + j] = val;
                }
            }
        }
    }

    template <typename TN>
    void DenseMatrix<TN>::resize(int64 rows, int64 cols, const TN& def_value)
    {
        std::vector<std::vector<TN>> new_data(rows, std::vector<TN>(cols, def_value));

        // 
        int64 min_rows = std::min(rows_, rows);
        int64 min_cols = std::min(columns_, cols);

        for (int64 i = 0; i < min_rows; ++i) {
            for (int64 j = 0; j < min_cols; ++j) {
                new_data[i][j] = c_matrix_[i * row_stride_ + j];
            }
        }

        // c_matrix_ = std::move(new_data);
        for (int64 i = 0; i < rows; ++i)
        {
            for (int64 j = 0; j < cols; ++j)
            {
                c_matrix_[i * row_stride_ + j] = new_data[i][j];
            }
        }

        rows_ = rows;
        columns_ = cols;
    }

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(const SparseMatrix<TN>& other)
        : rows_(other.rows()), columns_(other.cols()), 
        row_stride_(other.cols()), size_(other.rows() * other.cols()), 
        memory_owner_(true), c_matrix_(nullptr)
    {
        // Allocate memory for the dense matrix
        _allocate(size_);

        // Initialize all elements to zero
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] = TN(0);
        }

        // Populate non-zero elements from the SparseMatrix using CRS format
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 k = other.crs_row_ptr_[i]; k < other.crs_row_ptr_[i + 1]; ++k) {
                int64 j = other.crs_col_idx_[k]; // Column index
                TN value = other.crs_values_[k]; // Value at (i, j)
                c_matrix_[i * row_stride_ + j] = value;
            }
        }
    }

    template<typename TN>
    inline DenseMatrix<TN> DenseMatrix<TN>::submatrix(const std::pair<int64, int64>& rows_start_end,
        const std::pair<int64, int64>& cols_start_end)
    {
        int64 start_row = rows_start_end.first;
        int64 end_row = rows_start_end.second;
        int64 start_col = cols_start_end.first;
        int64 end_col = cols_start_end.second;

        if (start_row < 0 || start_row > end_row || end_row > rows_ ||
        start_col < 0 || start_col > end_col || end_col > columns_) 
        {
            throw std::out_of_range("Submatrix indices out of bounds");
        }

        int64 new_rows = end_row - start_row;
        int64 new_cols = end_col - start_col;
        DenseMatrix<TN> result(new_rows, new_cols);

        for (int64 i = 0; i < new_rows; ++i) 
        {
            for (int64 j = 0; j < new_cols; ++j) 
            {
                result.c_matrix_[i * result.row_stride_ + j] = c_matrix_[(start_row + i) * row_stride_ + (start_col + j)];
            }
        }

        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::identity(int64 size)
    {
        if (size < 0) 
        {
            throw std::invalid_argument("Size must be non-negative for identity matrix");
        }
        DenseMatrix<TN> result(size, size, TN(0));
        for (int64 i = 0; i < size; ++i) {
            result(i, i) = TN(1);
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::diagonal(const std::vector<TN>& values)
    {
        int64 size = static_cast<int64>(values.size());
        DenseMatrix<TN> result(size, size, TN(0));
        for (int64 i = 0; i < size; ++i) 
        {
            result(i, i) = values[i];
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::T(bool inplace)
    {
        DenseMatrix<TN> result(columns_, rows_);
        for (int i = 0; i < rows_; ++i)
        {
            for (int j = 0; j < columns_; ++j)
            {
                result(j, i) = (*this)(i, j);
            }
        }

        if (inplace)
        {
            *this = std::move(result);
        }

        return result;
    }

    template <typename TN>
    std::vector<TN> DenseMatrix<TN>::v() const
    {
        std::vector<TN> result(size_);
        std::copy_n(c_matrix_, size_, result.begin());
        return result;
    }

    template <typename TN>
    std::vector<std::vector<TN>> DenseMatrix<TN>::dd() const
    {
        std::vector<std::vector<TN>> result;
        result.reserve(rows_);

        for (int64 i = 0; i < rows_; ++i)
        {
            std::vector<TN> row(columns_);

            for (int64 j = 0; j < columns_; ++j)
            {
                row[j] = (*this)(i, j);
            }
            result.push_back(std::move(row));
        }

        return result;
    }

    template <typename TN>
    void DenseMatrix<TN>::_allocate(int64 capacity)
    {
        if (capacity < 0)
            throw std::invalid_argument("Capacity must be non-negative");

        capacity_ = std::max(static_cast<int64>(capacity * 1.5), static_cast<int64>(1));
        if (memory_owner_ && c_matrix_) 
        {
            TN* new_matrix = new TN[capacity_];
            for (int i = 0; i < size_; ++i)
                new_matrix[i] = c_matrix_[i];

            delete[] c_matrix_;
            c_matrix_ = new_matrix;
        }
        else if (!c_matrix_)
            c_matrix_ = new TN[capacity_]; 

        memory_owner_ = true;
    };

    template <typename TN>
    void DenseMatrix<TN>::_allocate(int64 rows, int64 cols)
    {
        if (rows < 0 || cols < 0) 
            throw std::invalid_argument("Rows and columns must be non-negative");

        capacity_ = static_cast<int64>(rows * cols * 1.2);
        _allocate(capacity_);
    };

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix() :
                    capacity_(1), size_(1), rows_(1),
                    columns_(1), row_stride_(1),
                    memory_owner_(true), c_matrix_(nullptr) {};

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(const int64& rows, const int64& columns, const TN& default_value) :
                    rows_(rows), columns_(columns), row_stride_(columns), size_(rows * columns),
                    c_matrix_(nullptr)
    {
        _allocate(rows * columns);
        for (int64 i = 0; i < size_; ++i)
            c_matrix_[i] = default_value;
    }

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(const TN* c_matrix, const int64& rows, const int64& columns)
                : rows_(rows), columns_(columns), row_stride_(columns), size_(rows * columns),
                memory_owner_(false), c_matrix_(const_cast<TN*>(c_matrix)) 
    {
        if (rows < 0 || columns < 0) {
            throw std::invalid_argument("Rows and columns must be non-negative");
        }
        if (!c_matrix && size_ > 0) {
            throw std::invalid_argument("Null pointer provided for non-empty matrix");
        }
        capacity_ = size_; 
    }

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(const std::vector<TN> c_matrix, const int64& rows, const int64& columns)
                    : rows_(rows), columns_(columns), row_stride_(columns), size_(rows * columns),
                    memory_owner_(true), c_matrix_(nullptr) 
    {
        if (rows < 0 || columns < 0) {
            throw std::invalid_argument("Rows and columns must be non-negative");
        }
        if (static_cast<int64>(c_matrix.size()) < size_) {
            throw std::invalid_argument("Vector size is smaller than required matrix size");
        }
        _allocate(rows, columns);
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] = c_matrix[i];
        }
    }

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(const std::vector<std::vector<TN>>& matrix)
                : rows_(matrix.size()), columns_(matrix.empty() ? 0 : matrix[0].size()),
                row_stride_(columns_), size_(rows_ * columns_),
                memory_owner_(true), c_matrix_(nullptr) 
    {
        if (rows_ == 0 && columns_ != 0) {
            throw std::invalid_argument("Empty matrix with non-zero columns is invalid");
        }
        for (int64 i = 0; i < rows_; ++i) {
            if (static_cast<int64>(matrix[i].size()) != columns_) {
                throw std::invalid_argument("All rows must have the same number of columns");
            }
        }
        _allocate(rows_, columns_); 
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                c_matrix_[i * row_stride_ + j] = matrix[i][j];
            }
        }
    }

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(const DenseMatrix<TN>& other)
                    : rows_(other.rows_), columns_(other.columns_), capacity_(other.capacity_),
                    size_(other.size_), row_stride_(other.row_stride_),
                    memory_owner_(true), c_matrix_(nullptr) 
    {
        _allocate(size_); 
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] = other.c_matrix_[i];
        }
    }

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(DenseMatrix<TN>&& other) noexcept
                    : rows_(other.rows_), columns_(other.columns_), capacity_(other.capacity_),
                    size_(other.size_), row_stride_(other.row_stride_),
                    memory_owner_(other.memory_owner_), c_matrix_(other.c_matrix_) {
        // Transfer ownership of memory and reset the source object
        other.c_matrix_ = nullptr;
        other.rows_ = 0;
        other.columns_ = 0;
        other.size_ = 0;
        other.capacity_ = 0;
        other.memory_owner_ = false;
    }

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(std::initializer_list<std::initializer_list<TN>> matrix)
                    : rows_(matrix.size()), columns_(matrix.begin()->size()), row_stride_(columns_),
                    size_(rows_ * columns_), memory_owner_(true), c_matrix_(nullptr) {
        for (const auto& row : matrix) {
            if (row.size() != columns_) {
                throw std::invalid_argument("All rows in the initializer list must have the same number of columns");
            }
        }
        
        _allocate(rows_, columns_);
        auto it = matrix.begin();
        for (int64 i = 0; i < rows_; ++i, ++it) {
            std::copy(it->begin(), it->end(), c_matrix_ + i * row_stride_);
        }
    }

    template <typename TN>
    DenseMatrix<TN>::DenseMatrix(std::initializer_list<TN> matrix, int rows, int columns)
                    : rows_(rows), columns_(columns), row_stride_(columns), size_(rows * columns),
                    memory_owner_(true), c_matrix_(nullptr) {

        if (rows < 0 || columns < 0) {
            throw std::invalid_argument("Rows and columns must be non-negative");
        }
        if (static_cast<int64>(matrix.size()) != size_) {
            throw std::invalid_argument("Initializer list size does not match the specified matrix size");
        }
        
        _allocate(rows, columns);
        std::copy(matrix.begin(), matrix.end(), c_matrix_);
    }

    template <typename TN>
    DenseMatrix<TN>::~DenseMatrix()
    {
        if (memory_owner_)
            delete[] c_matrix_;
    }

    template <typename TN>
    TN& DenseMatrix<TN>::operator()(int64 i, int64 j) 
    {
        if (i < 0 || i >= rows_ || j < 0 || j >= columns_) 
        {
            throw std::out_of_range("Matrix indices out of bounds");
        }
        return c_matrix_[i * row_stride_ + j];
    }

    template <typename TN>
    const TN& DenseMatrix<TN>::operator()(int64 i, int64 j) const 
    {
        if (i < 0 || i >= rows_ || j < 0 || j >= columns_) 
        {
            throw std::out_of_range("Matrix indices out of bounds");
        }
        return c_matrix_[i * row_stride_ + j];
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator+(const TN& scalar) const 
    {
        DenseMatrix<TN> result(*this); 
        for (int64 i = 0; i < size_; ++i) 
        {
            result.c_matrix_[i] += scalar; 
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator+(const TN*& c_matrix) const {
        DenseMatrix<TN> result(rows_, columns_); 
        for (int64 i = 0; i < size_; ++i) {
            result.c_matrix_[i] = c_matrix_[i] + c_matrix[i]; 
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator+(const TN** matrix) const {
        DenseMatrix<TN> result(rows_, columns_); 
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                result.c_matrix_[i * row_stride_ + j] = c_matrix_[i * row_stride_ + j] + matrix[i][j];
            }
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator+(const std::vector<TN> c_matrix) const {
        if (static_cast<int64>(c_matrix.size()) != size_) {
            throw std::invalid_argument("Vector size does not match matrix size");
        }
        DenseMatrix<TN> result(rows_, columns_); 
        for (int64 i = 0; i < size_; ++i) {
            result.c_matrix_[i] = c_matrix_[i] + c_matrix[i]; 
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator+(const std::vector<std::vector<TN>> c_matrix) const {
        if (static_cast<int64>(c_matrix.size()) != rows_ || 
            (rows_ > 0 && static_cast<int64>(c_matrix[0].size()) != columns_)) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        DenseMatrix<TN> result(rows_, columns_);
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                result.c_matrix_[i * row_stride_ + j] = c_matrix_[i * row_stride_ + j] + c_matrix[i][j];
            }
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator+(const DenseMatrix<TN>& other) const 
    {
        if (rows_ != other.rows_ || columns_ != other.columns_) 
        {
            throw std::invalid_argument("Matrix dimensions do not match for addition");
        }
        DenseMatrix<TN> result(rows_, columns_); 
        for (int64 i = 0; i < size_; ++i) 
        {
            result.c_matrix_[i] = c_matrix_[i] + other.c_matrix_[i];
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator-(const TN& scalar) const 
    {
        DenseMatrix<TN> result(*this); 
        for (int64 i = 0; i < size_; ++i) 
        {
            result.c_matrix_[i] -= scalar; 
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator-(const TN*& c_matrix) const 
    {
        DenseMatrix<TN> result(rows_, columns_); 
        for (int64 i = 0; i < size_; ++i) 
        {
            result.c_matrix_[i] = c_matrix_[i] - c_matrix[i]; 
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator-(const TN** matrix) const 
    {
        DenseMatrix<TN> result(rows_, columns_); 
        for (int64 i = 0; i < rows_; ++i) 
        {
            for (int64 j = 0; j < columns_; ++j) 
            {
                result.c_matrix_[i * row_stride_ + j] = c_matrix_[i * row_stride_ + j] - matrix[i][j];
            }
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator-(const std::vector<TN> c_matrix) const 
    {
        if (static_cast<int64>(c_matrix.size()) != size_) 
        {
            throw std::invalid_argument("Vector size does not match matrix size");
        }
        DenseMatrix<TN> result(rows_, columns_); 
        for (int64 i = 0; i < size_; ++i) 
        {
            result.c_matrix_[i] = c_matrix_[i] - c_matrix[i]; 
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator-(const std::vector<std::vector<TN>> c_matrix) const 
    {
        if (static_cast<int64>(c_matrix.size()) != rows_ || 
            (rows_ > 0 && static_cast<int64>(c_matrix[0].size()) != columns_)) 
        {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        DenseMatrix<TN> result(rows_, columns_);
        for (int64 i = 0; i < rows_; ++i) 
        {
            for (int64 j = 0; j < columns_; ++j) 
            {
                result.c_matrix_[i * row_stride_ + j] = c_matrix_[i * row_stride_ + j] - c_matrix[i][j];
            }
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator-(const DenseMatrix<TN>& other) const 
    {
        if (rows_ != other.rows_ || columns_ != other.columns_) 
        {
            throw std::invalid_argument("Matrix dimensions do not match for addition");
        }
        DenseMatrix<TN> result(rows_, columns_); 
        for (int64 i = 0; i < size_; ++i) 
        {
            result.c_matrix_[i] = c_matrix_[i] - other.c_matrix_[i];
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator*(const TN& scalar) const
    {
        DenseMatrix<TN> result(*this);
        for (int64 i = 0; i < size_; ++i) {
            result.c_matrix_[i] *= scalar;
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator*(const TN*& c_matrix) const
    {
        DenseMatrix<TN> result(rows_, columns_);
        for (int64 i = 0; i < size_; ++i) {
            result.c_matrix_[i] = c_matrix_[i] * c_matrix[i];
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator*(const TN** matrix) const
    {
        DenseMatrix<TN> result(rows_, columns_);
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                result.c_matrix_[i * row_stride_ + j] = c_matrix_[i * row_stride_ + j] * matrix[i][j];
            }
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator*(const std::vector<TN> c_matrix) const
    {
        if (static_cast<int64>(c_matrix.size()) != size_) {
            throw std::invalid_argument("Vector size does not match matrix size");
        }
        DenseMatrix<TN> result(rows_, columns_);
        for (int64 i = 0; i < size_; ++i) {
            result.c_matrix_[i] = c_matrix_[i] * c_matrix[i];
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator*(const std::vector<std::vector<TN>> c_matrix) const
    {
        if (static_cast<int64>(c_matrix.size()) != rows_ || 
            (rows_ > 0 && static_cast<int64>(c_matrix[0].size()) != columns_)) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        DenseMatrix<TN> result(rows_, columns_);
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                result.c_matrix_[i * row_stride_ + j] = c_matrix_[i * row_stride_ + j] * c_matrix[i][j];
            }
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator*(const DenseMatrix<TN>& other) const
    {
        if (columns_ != other.rows_) {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        DenseMatrix<TN> result(rows_, other.columns_);
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < other.columns_; ++j) {
                TN sum = 0;
                for (int64 k = 0; k < columns_; ++k) {
                    sum += c_matrix_[i * row_stride_ + k] * other.c_matrix_[k * other.row_stride_ + j];
                }
                result.c_matrix_[i * result.row_stride_ + j] = sum;
            }
        }
        return result;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator+=(const TN& scalar)
    {
        for (int64 i = 0; i < size_; ++i) 
        {
            c_matrix_[i] += scalar;
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator+=(const TN*& c_matrix)
    {
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] += c_matrix[i];
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator+=(const TN** matrix)
    {
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                c_matrix_[i * row_stride_ + j] += matrix[i][j];
            }
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator+=(const std::vector<TN> c_matrix)
    {
        if (static_cast<int64>(c_matrix.size()) != size_) {
            throw std::invalid_argument("Vector size does not match matrix size");
        }
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] += c_matrix[i];
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator+=(const std::vector<std::vector<TN>> c_matrix)
    {
        if (static_cast<int64>(c_matrix.size()) != rows_ || 
            (rows_ > 0 && static_cast<int64>(c_matrix[0].size()) != columns_)) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                c_matrix_[i * row_stride_ + j] += c_matrix[i][j];
            }
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator+=(const DenseMatrix<TN>& other)
    {
        if (rows_ != other.rows_ || columns_ != other.columns_) {
            throw std::invalid_argument("Matrix dimensions do not match for addition");
        }
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] += other.c_matrix_[i];
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator-=(const TN& scalar)
    {
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] -= scalar;
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator-=(const TN*& c_matrix)
    {
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] -= c_matrix[i];
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator-=(const TN** matrix)
    {
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                c_matrix_[i * row_stride_ + j] -= matrix[i][j];
            }
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator-=(const std::vector<TN> c_matrix)
    {
        if (static_cast<int64>(c_matrix.size()) != size_) {
            throw std::invalid_argument("Vector size does not match matrix size");
        }
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] -= c_matrix[i];
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator-=(const std::vector<std::vector<TN>> c_matrix)
    {
        if (static_cast<int64>(c_matrix.size()) != rows_ || 
            (rows_ > 0 && static_cast<int64>(c_matrix[0].size()) != columns_)) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                c_matrix_[i * row_stride_ + j] -= c_matrix[i][j];
            }
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator-=(const DenseMatrix<TN>& other)
    {
        if (rows_ != other.rows_ || columns_ != other.columns_) {
            throw std::invalid_argument("Matrix dimensions do not match for subtraction");
        }
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] -= other.c_matrix_[i];
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator*=(const TN& scalar)
    {
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] *= scalar;
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator*=(const TN*& c_matrix)
    {
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] *= c_matrix[i];
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator*=(const TN** matrix)
    {
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                c_matrix_[i * row_stride_ + j] *= matrix[i][j];
            }
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator*=(const std::vector<TN> c_matrix)
    {
        if (static_cast<int64>(c_matrix.size()) != size_) {
            throw std::invalid_argument("Vector size does not match matrix size");
        }
        for (int64 i = 0; i < size_; ++i) {
            c_matrix_[i] *= c_matrix[i];
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator*=(const std::vector<std::vector<TN>> c_matrix)
    {
        if (static_cast<int64>(c_matrix.size()) != rows_ || 
            (rows_ > 0 && static_cast<int64>(c_matrix[0].size()) != columns_)) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < columns_; ++j) {
                c_matrix_[i * row_stride_ + j] *= c_matrix[i][j];
            }
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator*=(const DenseMatrix<TN>& other)
    {
        if (columns_ != other.rows_) {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        DenseMatrix<TN> temp(rows_, other.columns_);
        for (int64 i = 0; i < rows_; ++i) {
            for (int64 j = 0; j < other.columns_; ++j) {
                TN sum = 0;
                for (int64 k = 0; k < columns_; ++k) {
                    sum += c_matrix_[i * row_stride_ + k] * other.c_matrix_[k * other.row_stride_ + j];
                }
                temp.c_matrix_[i * temp.row_stride_ + j] = sum;
            }
        }
        *this = std::move(temp);
        return *this;
    }

    template <typename TN>
    bool DenseMatrix<TN>::operator==(const DenseMatrix<TN>& other) const
    {
        if (rows_ != other.rows_ || columns_ != other.columns_) 
        {
            return false; 
        }
        for (int64 i = 0; i < size_; ++i) {
            if (c_matrix_[i] != other.c_matrix_[i]) 
            {
                return false; 
            }
        }
        return true; 
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator=(const DenseMatrix<TN>& other)
    {
        if (this != &other) { 
            if (memory_owner_ && c_matrix_) {
                delete[] c_matrix_; 
            }
            rows_ = other.rows_;
            columns_ = other.columns_;
            capacity_ = other.capacity_;
            size_ = other.size_;
            row_stride_ = other.row_stride_;
            memory_owner_ = true;
            c_matrix_ = new TN[capacity_];
            for (int64 i = 0; i < size_; ++i) {
                c_matrix_[i] = other.c_matrix_[i]; 
            }
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN>& DenseMatrix<TN>::operator=(DenseMatrix<TN>&& other) noexcept
    {
        if (this != &other) { 
            if (memory_owner_ && c_matrix_) {
                delete[] c_matrix_; 
            }
            rows_ = other.rows_;
            columns_ = other.columns_;
            capacity_ = other.capacity_;
            size_ = other.size_;
            row_stride_ = other.row_stride_;
            memory_owner_ = other.memory_owner_;
            c_matrix_ = other.c_matrix_; 
            
            other.c_matrix_ = nullptr;
            other.rows_ = 0;
            other.columns_ = 0;
            other.size_ = 0;
            other.capacity_ = 0;
            other.memory_owner_ = false;
        }
        return *this;
    }

    template <typename TN>
    DenseMatrix<TN> DenseMatrix<TN>::operator-() const
    {
        DenseMatrix<TN> result(*this);
        for (int64 i = 0; i < size_; ++i) {
            result.c_matrix_[i] = -result.c_matrix_[i];
        }
        return result;
    }
}