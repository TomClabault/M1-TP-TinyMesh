#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
public:
    Matrix() = delete;

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
     * matrix according to the given rotation angle parameters
     *
     * \param rX The rotation angle around the X axis
     * \param rY The rotation angle around the Y axis
     * \param rZ The rotation angle around the Z axis
     */
    void setRotation(double rX, double rY, double rZ);

    void setCoefficient(int y, int x, int coeff);

    int Rows() const;
    int Columns() const;

    Matrix Inverse();
    Matrix Tranpose();

    static Matrix RotationX(double xAngle);
    static Matrix RotationY(double yAngle);
    static Matrix RotationZ(double ZAngle);

private:
    const int MATRICE_HOMOTHETIE = 1;
    const int MATRICE_ROTATION = 2;

    int _matrixType;//This parameter is used when calculating the inverse / tranpose of
    //a matrix. It allows for a simple implementation (not a general one) of the
    //inverse and tranpose methods.

    int _rows, _columns;

    double** _coefficients;
};

#endif // MATRIX_H
