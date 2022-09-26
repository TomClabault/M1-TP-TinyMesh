#include "matrix.h"

#include <cmath>

Matrix::Matrix(double** values, int m, int n) : Matrix(m, n)
{
    for(int y = 0; y < m; y++)
        for(int x = 0; x < n; x++)
            _coefficients[y][x] = values[y][x];
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
}


void Matrix::setScaling(double sX, double sY, double sZ)
{
    _coefficients[0][0] = sX;
    _coefficients[1][1] = sY;
    _coefficients[2][2] = sZ;

    _matrixType = Matrix::MATRICE_HOMOTHETIE;
}

Matrix Matrix::RotationX(double xAngle)
{
    Matrix mat(3, 3);

    mat.setCoefficient(0, 0, 1);
    mat.setCoefficient(1, 1, std::cos(xAngle));
    mat.setCoefficient(1, 2, -std::sin(xAngle));
    mat.setCoefficient(2, 1, std::sin(xAngle));
    mat.setCoefficient(2, 2, std::cos(xAngle));

    mat._matr
    return mat;
}

static Matrix RotationY(double yAngle)
{
    Matrix mat(3, 3);

    mat.setCoefficient(1, 1, 1);
    mat.setCoefficient(0, 0, std::cos(yAngle));
    mat.setCoefficient(0, 2, std::sin(yAngle));
    mat.setCoefficient(2, 0, -std::sin(yAngle));
    mat.setCoefficient(2, 2, std::cos(yAngle));

    return mat;
}

static Matrix RotationZ(double ZAngle)
{
    Matrix mat(3, 3);

    mat.setCoefficient(2, 2, 1);
    mat.setCoefficient(0, 0, std::cos(ZAngle));
    mat.setCoefficient(0, 1, -std::sin(ZAngle));
    mat.setCoefficient(1, 0, std::sin(ZAngle));
    mat.setCoefficient(1, 1, std::cos(ZAngle));

    return mat;
}

void Matrix::setRotation(double rX, double rY, double rZ)
{
    Matrix matRX = Matrix::RotationX(rX);
    Matrix matRY = Matrix::RotationX(rY);
    Matrix matRZ = Matrix::RotationX(rZ);

    _coefficients[0][0] = sX;
    _coefficients[1][1] = sY;
    _coefficients[2][2] = sZ;

    _matrixType = Matrix::MATRICE_ROTATION;
}

void Matrix::setCoefficient(int y, int x, int coeff)
{
    _coefficients[y][x] = coeff;
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

}

Matrix Matrix::Tranpose()
{

}
