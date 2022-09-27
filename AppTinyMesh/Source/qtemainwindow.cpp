#include "qte.h"
#include "implicits.h"

#include "icosphere.h"

#include "matrix.h"

MainWindow::MainWindow()
{
//    /*!
//     * \brief Fills the coefficient of the matrix to
//     * have this instance of Matrix behave like a rotation
//     * matrix according to the given rotation angle parameters.
//     * The rotations are applied in the following order:
//     * Z -> Y -> X
//     *
//     * \param rX The rotation angle around the X axis in radians
//     * \param rY The rotation angle around the Y axis in radians
//     * \param rZ The rotation angle around the Z axis in radians
//     */
//    void setRotation(double rX, double rY, double rZ);

//    int Rows() const;
//    int Columns() const;

//    Matrix Inverse();
//    Matrix Tranpose();

//    double* operator()(int y) const;
//    double& operator()(int y, int x) const;

//    Matrix& operator=(const Matrix& A);

//    Matrix& operator+=(const Matrix& B);
//    Matrix& operator-=(const Matrix& B);
//    Matrix& operator*=(const Matrix& B);
//    Matrix& operator*=(double n);
//    Matrix& operator/=(const Matrix& B);
//    Matrix& operator/=(double n);

//    friend Matrix operator+(const Matrix& A, const Matrix& B);
//    friend Matrix operator-(const Matrix& A, const Matrix& B);
//    friend Matrix operator*(const Matrix& A, const Matrix& B);
//    friend Matrix operator*(const Matrix& A, double);
//    friend Matrix operator/(const Matrix& A, double);

//    static Matrix RotationX(double xAngle)
//    {
//        Matrix mat(3, 3);

//        mat(0, 0) = 1;
//        mat(1, 1) = std::cos(xAngle);
//        mat(1, 2) = -std::sin(xAngle);
//        mat(2, 1) = std::sin(xAngle);
//        mat(2, 2) = std::cos(xAngle);

//        return mat;
//    }

//    static Matrix RotationY(double yAngle)
//    {
//        Matrix mat(3, 3);

//        mat(1, 1) = 1;
//        mat(0, 0) = std::cos(yAngle);
//        mat(0, 2) = std::sin(yAngle);
//        mat(2, 0) = -std::sin(yAngle);
//        mat(2, 2) = std::cos(yAngle);

//        return mat;
//    }

//    static Matrix RotationZ(double ZAngle)
//    {
//        Matrix mat(3, 3);

//        mat(2, 2) = 1;
//        mat(0, 0) = std::cos(ZAngle);
//        mat(0, 1) = -std::sin(ZAngle);
//        mat(1, 0) = std::sin(ZAngle);
//        mat(1, 1) = std::cos(ZAngle);

//        return mat;
//    }

    double** arrayMat = new double*[5];
    for(int i = 0; i < 5; i++)
    {
        arrayMat[i] = new double[5];

        for(int j = 0; j < 5; j++)
            arrayMat[i][j] = i * 5 + j;
    }


    Matrix A(5, 5);
    A.setRotation(1, 2, 3);

    std::cout << A;

    std::exit(10);

    // Chargement de l'interface
	uiw.setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
	uiw.widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}

MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
	connect(uiw.boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw.sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw.icosphereMesh, SIGNAL(clicked()), this, SLOT(IcosphereMeshExample()));
	connect(uiw.resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
	connect(uiw.wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
	connect(uiw.radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
	connect(uiw.radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

void MainWindow::BoxMeshExample()
{
	Mesh boxMesh = Mesh(Box(1.0));

	std::vector<Color> cols;
	cols.resize(boxMesh.Vertexes());
	for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

	meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
	UpdateGeometry();
}

void MainWindow::SphereImplicitExample()
{
  AnalyticScalarField implicit;

  Mesh implicitMesh;
  implicit.Polygonize(31, implicitMesh, Box(2.0));

  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (int i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::IcosphereMeshExample()
{
    Mesh icosphereMesh = Mesh(Icosphere());

    std::vector<Color> colors;
    colors.resize(icosphereMesh.Vertexes());
    for (int i = 0; i < colors.size(); i++)
        colors[i] = Color(1.0, 0.5, 0.5);

    meshColor = MeshColor(icosphereMesh, colors, icosphereMesh.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh", meshColor);

	uiw.lineEdit->setText(QString::number(meshColor.Vertexes()));
	uiw.lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
	meshWidget->UseWireframeGlobal(uiw.wireframe->isChecked());

	if (uiw.radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
