#include "matrix.h"

#include <cmath>
#include <stdexcept>

Matrix::Matrix(const Matrix& A) : Matrix(A._coefficients, A._rows, A._columns)
{
    _matrixType = A._matrixType;
}

Matrix::Matrix(Matrix&& A)
{
    _coefficients = A._coefficients;

    A._coefficients = nullptr;

    _rows = A._rows;
    _columns = A._columns;

    _matrixType = A._matrixType;
}

Matrix::Matrix(double** values, int m, int n)
{
    _coefficients = values;

    _rows = m;
    _columns = n;

    _matrixType = -1;
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

void Matrix::setScaling(double sX, double sY, double sZ)
{
    _matrixType = Matrix::MATRICE_HOMOTHETIE;

    if(_rows != 3 || _columns != 3)//Resizing the matrix
        _reallocMem(3, 3);

    _coefficients[0][0] = sX;
    _coefficients[1][1] = sY;
    _coefficients[2][2] = sZ;
}

void Matrix::setRotation(double rX, double rY, double rZ)
{
    _matrixType = Matrix::MATRICE_ROTATION;

    if(_rows != 3 || _columns != 3)//Resizing the matrix
        _reallocMem(3, 3);

    (*this) = Matrix::RotationX(rX) * Matrix::RotationY(rY) * Matrix::RotationZ(rZ);
}

int Matrix::Rows() const
{
    return _rows;
}

int Matrix::Columns() const
{
    return _columns;
}

Matrix Matrix::Inverse()
{
    return Matrix(3, 3);
}

Matrix Matrix::Tranpose()
{
    return Matrix(3, 3);
}

double* Matrix::operator()(int y) const
{
    return this->_coefficients[y];
}

double& Matrix::operator()(int y, int x) const
{
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
