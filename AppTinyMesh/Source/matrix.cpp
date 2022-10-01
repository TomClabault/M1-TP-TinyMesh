#include "matrix.h"

#include <cmath>
#include <stdexcept>

Matrix::Matrix(const Matrix& A)
{
    _coefficients = new double*[A._rows];
    for(int i = 0; i < A._rows; i++)
    {
        _coefficients[i] = new double[A._columns];
        for(int j = 0; j < A._columns; j++)
            _coefficients[i][j] = A(i, j);
    }

    _rows = A._rows;
    _columns = A._columns;

    _matrixType = A._matrixType;
    _isTransposed = A._isTransposed;
}

Matrix::Matrix(Matrix&& A)
{
    _coefficients = A._coefficients;
    A._coefficients = nullptr;

    _rows = A._rows;
    _columns = A._columns;

    _matrixType = A._matrixType;
    _isTransposed = A._isTransposed;
}

Matrix::Matrix(double** values, int m, int n)
{
    _coefficients = values;

    _rows = m;
    _columns = n;
}

Matrix::Matrix(int m, int n)
{
    _coefficients = new double*[m];
    for(int i = 0; i < m; i++)
    {
        _coefficients[i] = new double[n];
        for(int j = 0; j < n; j++)
            _coefficients[i][j] = 0;
    }

    _rows = m;
    _columns = n;

    _matrixType = -1;
}

Matrix::~Matrix()
{
    _freeMem();
}

int Matrix::Rows() const
{
    return _rows;
}

int Matrix::Columns() const
{
    return _columns;
}

Matrix& Matrix::Transpose()
{
    if(_matrixType != HOMOTHETY_MATRIX && _matrixType != ROTATION_MATRIX)
        throw std::invalid_argument("Inverse only implemented for rotation and homothety matrices.");
    else
        this->_isTransposed = !this->_isTransposed;//This will allow to just swap the
        //indexes when accessing an element with the () operator

    return *this;
}

Matrix& Matrix::Inverse()
{
    if(_matrixType == ROTATION_MATRIX)//The inverse of a rotation matrix is its transposition
        //because rotation matrices are orthogonal
        return this->Transpose();
    else if(_matrixType == HOMOTHETY_MATRIX)
    {
        (*this)(0, 0) = 1 / (*this)(0, 0);
        (*this)(1, 1) = 1 / (*this)(1, 1);
        (*this)(2, 2) = 1 / (*this)(2, 2);

        return *this;
    }
    else
        throw std::invalid_argument("Inverse only implemented for rotation and homothety matrices.");
}

double& Matrix::operator()(int y, int x) const
{
    if(_isTransposed)
        std::swap(y, x);

    return this->_coefficients[y][x];
}

Matrix& Matrix::operator=(const Matrix& A)
{
    if(this->_rows != A._rows || this->_columns != A._columns)
        throw std::invalid_argument("The matrices must be the same size");

    for(int y = 0; y < this->_rows; y++)
        for(int x = 0; x < this->_rows; x++)
            (*this)(y, x) = A(y, x);

    return (*this);
}

Matrix& Matrix::operator=(Matrix&& A)
{
    this->_coefficients = A._coefficients;
    this->_rows = A._rows;
    this->_columns = A._columns;
    this->_matrixType = A._matrixType;

    A._coefficients = nullptr;

    return (*this);
}

Matrix operator+(const Matrix& A, const Matrix& B)
{
    if(A._rows != B._rows || A._columns != B._columns)
        throw std::invalid_argument("The matrices must be the same size");

    Matrix result(A._rows, A._columns);
    for(int y = 0; y < A._rows; y++)
        for(int x = 0; x < A._columns; x++)
            result(y, x) = A(y, x) + B(y, x);

    return result;
}

Matrix& Matrix::operator+=(const Matrix& B)
{
    return (*this) = ((*this) + B);
}


Matrix operator-(const Matrix& A, const Matrix& B)
{
    if(A._rows != B._rows || A._columns != B._columns)
        throw std::invalid_argument("The matrices must be the same size");

    Matrix result(A._rows, A._columns);
    for(int y = 0; y < A._rows; y++)
        for(int x = 0; x < A._columns; x++)
            result(y, x) = A(y, x) - B(y, x);

    return result;
}

Matrix& Matrix::operator-=(const Matrix& B)
{
    return (*this) = ((*this) - B);
}

Matrix operator*(const Matrix& A, const Matrix& B)
{
    if(A._columns != B._rows)
        throw std::invalid_argument("The matrices are not compatible");

    Matrix result(A._rows, B._columns);

    for(int y = 0; y < A._rows; y++)
        for(int x = 0; x < B._columns; x++)
            for(int i = 0; i < A._columns; i++)
                result(y, x) += A(y, i) * B(i, x);

    if(A._matrixType == B._matrixType)
        result._matrixType = A._matrixType;

    return result;
}

Matrix& Matrix::operator*=(const Matrix& B)
{
    Matrix temp(B);

    return (*this) = ((*this) * temp);
}

Matrix& Matrix::operator*=(double n)
{
    return (*this) = ((*this) * n);
}

Matrix operator*(const Matrix& A, double n)
{
    Matrix result(A);

    for(int y = 0; y < A._rows; y++)
        for(int x = 0; x < A._columns; x++)
            A(y, x) *= n;

    return result;
}

Matrix operator/(const Matrix& A, double n)
{
    Matrix result(A);

    for(int y = 0; y < A._rows; y++)
        for(int x = 0; x < A._columns; x++)
            A(y, x) /= n;

    return result;
}

Matrix& Matrix::operator/=(double n)
{
    return (*this) = ((*this) / n);
}

std::ostream& operator << (std::ostream& os, const Matrix& A)
{
    for(int y = 0; y < A._rows; y++)
    {
        for(int x = 0; x < A._columns; x++)
            os << A(y, x) << " ";

        os << std::endl;
    }

    return os;
}
