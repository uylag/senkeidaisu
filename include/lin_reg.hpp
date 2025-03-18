#pragma once
#include <senkeidaisu.hpp>
#include <iostream>
#include <initializer_list>
using namespace senkeidaisu;

template <typename TN>
class LinearRegression
{
public:
    LinearRegression(): X_(), y_(), theta_(), learning_rate_(TN(0)) {
        std::cout << "[LinearRegression] default constructor called.";
    };

    LinearRegression(const Matrix<TN>& X, const Matrix<TN>& y, const TN& learning_rate = 0.01):
                     X_(X), y_(y), theta_(Matrix<TN>(X.cols(), 1, TN(0))), learning_rate_(learning_rate) 
                     {
                        std::cout << "[LinearRegression] Initialized with:\n";
                        std::cout << "  X shape: " << X_.rows() << "x" << X_.cols() << "\n";
                        std::cout << "  y shape: " << y_.rows() << "x" << y_.cols() << "\n";
                        std::cout << "  Learning rate: " << learning_rate_ << "\n";
                     };  

    Matrix<TN> predict() 
    { 
        std::cout << "[predict] Called\n";
        std::cout << "  X shape: " << X_.rows() << "x" << X_.cols() << "\n";
        std::cout << "  theta shape: " << theta_.rows() << "x" << theta_.cols() << "\n";

        Matrix<TN> result = X_ * theta_;

        std::cout << "  Predicted output shape: " << result.rows() << "x" << result.cols() << "\n";
        return result; 
    };

    Matrix<TN> predict(std::initializer_list<TN> list)
    {
        std::cout << "[predict] with init list called.\n";

        if (theta_.rows() == 0 || theta_.cols() == 0) {
            throw std::runtime_error("[predict] Model is not fitted yet. Call fit() first.");
        }

        Matrix<TN> X_input(1, list.size());
        int col = 0;

        for (const TN& value: list)
        {
            X_input(0, col++) = value;
        }

        std::cout << "  Input shape: " << X_input.rows() << "x" << X_input.cols() << "\n";

        if (X_input.cols() != X_.cols()) {
            throw std::invalid_argument("[predict] Input feature count does not match trained model.");
        }

        Matrix<TN> result = X_input * theta_;
        std::cout << "Output shape: " << result.rows() << "x" << result.cols() << "\n";

        return result;
    }

    void fit(const Matrix<TN>& X, const Matrix<TN>& y)
    {
        std::cout << "[fit] Called\n";
        std::cout << "  X shape: " << X.rows() << "x" << X.cols() << "\n";
        std::cout << "  y shape: " << y.rows() << "x" << y.cols() << "\n";

        if (y.rows() == 1) {
            y_ = y.T(false); // Convert to column vector
        } else {
            y_ = y;
        }

        theta_ = grad_desc(100);
        std::cout << "  Updated theta shape: " << theta_.rows() << "x" << theta_.cols() << "\n";
    }

    TN MSE()
    {
        std::cout << "[MSE] called\n";

        int m = X_.rows();
        Matrix<TN> pred = predict();

        std::cout << "Predicted shape: " << pred.rows() << "x" << pred.cols() << "\n";
        std::cout << "Target shape: " << y_.rows() << "x" << y_.cols() << "\n";

        if (pred.rows() != y_.rows() || pred.cols() != y_.cols()) 
        {
            throw std::invalid_argument("[MSE] Dimension mismatch: predicted and target sizes are not equal.");
        }

        Matrix<TN> error = pred - y_;

        std::cout << "  Error matrix shape: " << error.rows() << "x" << error.cols() << "\n";

        TN cost = (error.T(false) * error)(0, 0) / (2 * m);

        std::cout << "  Computed MSE: " << cost << "\n";
        return cost;
    };

    Matrix<TN> grad_desc(int64 iterations) 
    {
        std::cout << "[grad_desc] Called\n";
        std::cout << "  X shape: " << X_.rows() << "x" << X_.cols() << "\n";
        std::cout << "  y shape: " << y_.rows() << "x" << y_.cols() << "\n";
        std::cout << "  Iterations: " << iterations << "\n";

        if (X_.cols() != theta_.rows()) 
        {
            throw std::invalid_argument("[grad_desc] Dimension mismatch: X_.cols() != theta_.rows()");
        }

        Matrix<TN> result = X_.gradient_descent(X_, y_, learning_rate_, iterations);
        std::cout << "  Gradient descent complete. New theta shape: " << result.rows() << "x" << result.cols() << "\n";
        
        return result;
    }

    Matrix<TN> get_theta()
    {
        std::cout << "[get_theta] Called\n";
        std::cout << "  X shape: " << X_.rows() << "x" << X_.cols() << "\n";

        if (X_.rows() < X_.cols()) 
        {
            throw std::invalid_argument("[get_theta] X_.rows() < X_.cols(), cannot compute normal equation.");
        }

        Matrix<TN> equation = (X_.T(false) * X_).inv() * X_.T(false) * y_;
        
        std::cout << "  Computed theta shape: " << equation.rows() << "x" << equation.cols() << "\n";

        theta_ = equation;
        theta_ = grad_desc(100);

        std::cout << "  Final theta shape: " << theta_.rows() << "x" << theta_.cols() << "\n";

        return theta_;
    }



private:
    Matrix<TN> X_;
    Matrix<TN> y_;
    Matrix<TN> theta_;
    TN learning_rate_;
};