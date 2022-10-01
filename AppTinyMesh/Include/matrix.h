#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

#include "mathematics.h"


class Matrix
{
public:
    Matrix() = delete;

    Matrix(const Matrix& A);
    Matrix(Matrix&& A);

    /*!
     * \brief Instantiates a matric of size m*x with
     * all coefficients initialized at 0
     *
     * \param m Number of rows of the matrix
     * \param n Number of columns of the matrix
    */
    Matrix(int m, int n);

    /*!
     * \brief Instantiates a matrix of size m * n and
     * fills the coefficient of the matrix with the values
     * found in the array passed as parameter
     *
     * \param values The values that are going to be used
     * to initialize the newly created matrix
     * \param m The number of rows of the matrix
     * \param n The number of columns of the matrix
     */
    Matrix(double** values, int m, int n);

    ~Matrix();

    int Rows() const;
    int Columns() const;

    Matrix& Inverse();
    Matrix GetInverse() const;

    Matrix& Transpose();
    Matrix GetTranspose() const;

    double& operator()(int y, int x) const;

    Matrix& operator=(const Matrix& A);
    Matrix& operator=(Matrix&& A);

    Matrix& operator+=(const Matrix& B);
    Matrix& operator-=(const Matrix& B);
    Matrix& operator*=(const Matrix& B);
    Matrix& operator*=(double n);
    Matrix& operator/=(double n);

    friend Matrix operator+(const Matrix& A, const Matrix& B);
    friend Matrix operator-(const Matrix& A, const Matrix& B);
    friend Vector operator*(const Vector& A, const Matrix& B);
    friend Matrix operator*(const Matrix& A, const Matrix& B);
    friend Matrix operator*(const Matrix& A, double);
    friend Matrix operator/(const Matrix& A, double);

    static Matrix RotationX(double xAngle)
    {
        Matrix mat(3, 3);

        mat(0, 0) = 1;
        mat(1, 1) = std::cos(xAngle);
        mat(1, 2) = -std::sin(xAngle);
        mat(2, 1) = std::sin(xAngle);
        mat(2, 2) = std::cos(xAngle);

        mat._matrixType = ROTATION_MATRIX;

        return mat;
    }

    static Matrix RotationY(double yAngle)
    {
        Matrix mat(3, 3);

        mat(1, 1) = 1;
        mat(0, 0) = std::cos(yAngle);
        mat(0, 2) = std::sin(yAngle);
        mat(2, 0) = -std::sin(yAngle);
        mat(2, 2) = std::cos(yAngle);

        mat._matrixType = ROTATION_MATRIX;

        return mat;
    }

    static Matrix RotationZ(double ZAngle)
    {
        Matrix mat(3, 3);

        mat(2, 2) = 1;
        mat(0, 0) = std::cos(ZAngle);
        mat(0, 1) = -std::sin(ZAngle);
        mat(1, 0) = std::sin(ZAngle);
        mat(1, 1) = std::cos(ZAngle);

        mat._matrixType = ROTATION_MATRIX;

        return mat;
    }

    /*!
     * \brief Creates and returns a rotation matrix that applies
     * 3 rotations depending on the given angle parameters. The
     * rotation in the Z axis is first applied, then Y, then X.
     *
     * \param xAngle The angle for the rotation around the X axis.
     * Angle in radians
     * \param yAngle The angle for the rotation around the Y axis.
     * Angle in radians
     * \param zAngle The angle for the rotation around the Z axis.
     * Angle in radians
     *
     * \return The rotation matrix that applies the rotations in
     * the order specified in the brief of the function
     */
    static Matrix RotationZYX(double xAngle, double yAngle, double zAngle)
    {
        Matrix mat = RotationX(xAngle) * RotationY(yAngle) * RotationZ(zAngle);

        return mat;
    }

    /*!
     * \brief Overload for RotationZYX(angle, angle, angle)
     *
     * \param angle The angle around which to apply the rotation
     * for the Z, Y and X axis. Angle in radians.
     *
     * \return A rotation matrix that rotates around all the axis
     * by the given angle. Rotation around Z first, then Y, then X
     */
    static Matrix RotationZYX(double angle)
    {
        return RotationZYX(angle, angle, angle);
    }

    static Matrix Homothety(double xScale, double yScale, double zScale)
    {
        Matrix mat(3,3 );

        mat(0, 0) = xScale;
        mat(1, 1) = yScale;
        mat(2, 2) = zScale;

        mat._matrixType = HOMOTHETY_MATRIX;

        return mat;
    }

    static Matrix Homothety(Vector s)
    {
        return Homothety(s[0], s[1], s[2]);
    }

    friend std::ostream& operator << (std::ostream& os, const Matrix& A);

private:
    const static int HOMOTHETY_MATRIX = 1;
    const static int ROTATION_MATRIX = 2;

    bool _isTransposed = false;

    int _matrixType = -1;//This parameter is used when calculating the inverse / tranpose of
    //a matrix. It allows for a simple implementation (not a general one) of the
    //inverse and tranpose methods.

    int _rows, _columns;

    double** _coefficients = nullptr;

    void _freeMem()
    {
        //The memory has already been moved to another location
        //We don't have to free anything
        if(_coefficients == nullptr)
            return;

        for(int i = 0; i < this->_rows; i++)
            delete[] _coefficients[i];

        delete _coefficients;
    }
};

#endif // MATRIX_H
