#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

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

    /*!
     * \brief Fills the coefficient of the matrix to
     * have this instance of Matrix behave like a scaling
     * matrix according to the given scaling parameters
     *
     * \param sX The X axis scaling
     * \param sY The Y axis scaling
     * \param sZ The Z axis scaling
     */
    void setScaling(double sX, double sY, double sZ);

    /*!
     * \brief Fills the coefficient of the matrix to
     * have this instance of Matrix behave like a rotation
     * matrix according to the given rotation angle parameters.
     * The rotations are applied in the following order:
     * Z -> Y -> X
     *
     * \param rX The rotation angle around the X axis in radians
     * \param rY The rotation angle around the Y axis in radians
     * \param rZ The rotation angle around the Z axis in radians
     */
    void setRotation(double rX, double rY, double rZ);

    int Rows() const;
    int Columns() const;

    Matrix Inverse();
    Matrix Tranpose();

    double* operator()(int y) const;
    double& operator()(int y, int x) const;

    Matrix& operator=(const Matrix& A);

    Matrix& operator+=(const Matrix& B);
    Matrix& operator-=(const Matrix& B);
    Matrix& operator*=(const Matrix& B);
    Matrix& operator*=(double n);
    Matrix& operator/=(const Matrix& B);
    Matrix& operator/=(double n);

    friend Matrix operator+(const Matrix& A, const Matrix& B);
    friend Matrix operator-(const Matrix& A, const Matrix& B);
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

        return mat;
    }

    friend std::ostream& operator << (std::ostream& os, const Matrix& A);

private:
    const int MATRICE_HOMOTHETIE = 1;
    const int MATRICE_ROTATION = 2;

    int _matrixType;//This parameter is used when calculating the inverse / tranpose of
    //a matrix. It allows for a simple implementation (not a general one) of the
    //inverse and tranpose methods.

    int _rows, _columns;

    double** _coefficients = nullptr;

    void _freeMem()
    {
        for(int i = 0; i < this->_rows; i++)
            delete this->_coefficients[i];

        delete this->_coefficients;
    }

    void inline _reallocMem(int m, int n)
    {
        _freeMem();

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
};

#endif // MATRIX_H
