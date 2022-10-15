#include "qtemainwindow.h"
#include "simpleMeshes.h"
#include "implicits.h"
#include "ui_interface.h"

MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
	// Chargement de l'interface
    uiw->setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
    uiw->widget_GL->setLayout(GLlayout);

    uiw->toolboxGroupBox->setVisible(false);

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
    connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw->icosphereButton, SIGNAL(clicked()), this, SLOT(DisplayIcosphere()));
    connect(uiw->torusButton, SIGNAL(clicked()), this, SLOT(DisplayTorus()));
    connect(uiw->capsuleButton, SIGNAL(clicked()), this, SLOT(DisplayCapsule()));
    connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
    connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

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
    uiw->toolboxGroupBox->setVisible(false);

	Mesh boxMesh = Mesh(Box(1.0));

	std::vector<Color> cols;
	cols.resize(boxMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

	meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
	UpdateGeometry();
}

void MainWindow::SphereImplicitExample()
{
    uiw->toolboxGroupBox->setVisible(false);

    AnalyticScalarField implicit;

    Mesh implicitMesh;
    implicit.Polygonize(31, implicitMesh, Box(2.0));

    std::vector<Color> cols;
    cols.resize(implicitMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(0.8, 0.8, 0.8);

    meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::CreateIcosphereMesh(double radius, int subdivisions)
{
    Mesh icosphereMesh = Mesh(Icosphere(radius, subdivisions));

    std::vector<Color> cols;
    cols.resize(icosphereMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(icosphereMesh, cols, icosphereMesh.VertexIndexes());
}

void MainWindow::CreateTorusMesh(double innerRadius, double outerRadius, int ringCount, int ringsSubdivisions)
{
    //Mesh torusMesh = Mesh(Torus(innerRadius, outerRadius, ringCount, ringsSubdivisions));
    Mesh torusMesh = Mesh(Cylinder(1, 2, 1, 10));

    std::vector<Color> cols;
    cols.resize(torusMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(torusMesh, cols, torusMesh.VertexIndexes());
}

void MainWindow::CreateCapsuleMesh(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions)
{
    Mesh capsuleMesh = Mesh(Capsule(radius, cylinderHeight, cylinderHeightSubdivions, cylinderSubdivisions, sphereHeightSubdivisions));

    std::vector<Color> cols;
    cols.resize(capsuleMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(capsuleMesh, cols, capsuleMesh.VertexIndexes());
}

void MainWindow::SetupIcosphereToolbox()
{
    uiw->toolboxGroupBox->setVisible(true);

    delete toolboxWidget;//Deleting the previous widget
    toolboxWidget = new QWidget;
    icosphereToolbox.setupUi(toolboxWidget);

    QVBoxLayout vBoxLayout(uiw->toolboxGroupBox);
    vBoxLayout.addWidget(toolboxWidget);

    //Default settings for the icosphere at its creation
    icosphereToolbox.icosphereRadiusInput->setText("1.0");
    icosphereToolbox.icosphereSubdivisionInput->setText("1");

    connect(icosphereToolbox.icosphereRadiusInput, SIGNAL(returnPressed()), this, SLOT(UpdateIcosphere()));
    connect(icosphereToolbox.icosphereSubdivisionInput, SIGNAL(returnPressed()), this, SLOT(UpdateIcosphere()));
    connect(icosphereToolbox.applyIcosphereToolboxButton, SIGNAL(clicked()), this, SLOT(UpdateIcosphere()));
}

void MainWindow::SetupTorusToolbox()
{
    uiw->toolboxGroupBox->setVisible(true);

    delete toolboxWidget;//Deleting the previous widget
    toolboxWidget = new QWidget;
    torusToolbox.setupUi(toolboxWidget);

    QVBoxLayout vboxLayout(uiw->toolboxGroupBox);
    vboxLayout.addWidget(toolboxWidget);

    //Default settings for the icosphere at its creation
    torusToolbox.torusInRadiusInput->setText("0.375");
    torusToolbox.torusOutRadiusInput->setText("1.5");
    torusToolbox.torusRingCountInput->setText("20");
    torusToolbox.torusRingSubdivInput->setText("20");

    connect(torusToolbox.torusInRadiusInput, SIGNAL(returnPressed()), this, SLOT(UpdateTorus()));
    connect(torusToolbox.torusOutRadiusInput, SIGNAL(returnPressed()), this, SLOT(UpdateTorus()));
    connect(torusToolbox.torusRingCountInput, SIGNAL(returnPressed()), this, SLOT(UpdateTorus()));
    connect(torusToolbox.torusRingSubdivInput, SIGNAL(returnPressed()), this, SLOT(UpdateTorus()));
    connect(torusToolbox.applyTorusToolboxButton, SIGNAL(clicked()), this, SLOT(UpdateTorus()));
}

void MainWindow::SetupCapsuleToolbox()
{
    uiw->toolboxGroupBox->setVisible(true);

    delete toolboxWidget;//Deleting the previous widget
    toolboxWidget = new QWidget;
    capsuleToolbox.setupUi(toolboxWidget);

    QVBoxLayout vboxLayout(uiw->toolboxGroupBox);
    vboxLayout.addWidget(toolboxWidget);

    //Default settings for the icosphere at its creation
    capsuleToolbox.capsuleRadiusInput->setText("1.0");
    capsuleToolbox.capsuleCylinderHeightInput->setText("2.0");
    capsuleToolbox.capsuleCylinderHeightSubdivInput->setText("5");
    capsuleToolbox.capsuleCylinderSubdivInput->setText("10");
    capsuleToolbox.capsuleCapsSubdivInput->setText("10");

    connect(capsuleToolbox.capsuleRadiusInput, SIGNAL(returnPressed()), this, SLOT(UpdateCapsule()));
    connect(capsuleToolbox.capsuleCylinderHeightInput, SIGNAL(returnPressed()), this, SLOT(UpdateCapsule()));
    connect(capsuleToolbox.capsuleCylinderHeightSubdivInput, SIGNAL(returnPressed()), this, SLOT(UpdateCapsule()));
    connect(capsuleToolbox.capsuleCylinderSubdivInput, SIGNAL(returnPressed()), this, SLOT(UpdateCapsule()));
    connect(capsuleToolbox.capsuleCapsSubdivInput, SIGNAL(returnPressed()), this, SLOT(UpdateCapsule()));
    connect(capsuleToolbox.applyCapsuleToolboxButton, SIGNAL(clicked()), this, SLOT(UpdateCapsule()));
}

void MainWindow::DisplayIcosphere()
{
    CreateIcosphereMesh(1.0, 1);

    UpdateGeometry();

    SetupIcosphereToolbox();
}

void MainWindow::DisplayTorus()
{
    CreateTorusMesh(0.375, 1.5, 20, 20);

    UpdateGeometry();

    SetupTorusToolbox();
}

void MainWindow::DisplayCapsule()
{
    CreateCapsuleMesh(1.0, 2.0, 5, 10, 10);

    UpdateGeometry();

    SetupCapsuleToolbox();
}

void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh", meshColor);

    uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

    if (uiw->radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}

double getSafeDoubleFromInput(const QLineEdit* input)
{
    double value;

    QString inputText = input->text();
    try
    {
        value = std::stod(inputText.toStdString());
    }
    catch (std::invalid_argument e)
    {
        return -1;
    }
    catch (std::out_of_range e)
    {
        return -1;
    }

    return value;
}

int getSafeIntFromInput(const QLineEdit* input)
{
    int value;

    QString inputText = input->text();
    try
    {
        value = std::stoi(inputText.toStdString());
    }
    catch (std::invalid_argument e)
    {
        return -1;
    }
    catch (std::out_of_range e)
    {
        return -1;
    }

    return value;
}

double MainWindow::getIcosphereToolboxRadius()
{
    return getSafeDoubleFromInput(icosphereToolbox.icosphereRadiusInput);
}

int MainWindow::getIcosphereToolboxSubdiv()
{
    return getSafeIntFromInput(icosphereToolbox.icosphereSubdivisionInput);
}

double MainWindow::getTorusToolboxInRadius()
{
    return getSafeDoubleFromInput(torusToolbox.torusInRadiusInput);
}

double MainWindow::getTorusToolboxOutRadius()
{
    return getSafeDoubleFromInput(torusToolbox.torusOutRadiusInput);
}

int MainWindow::getTorusToolboxRingCount()
{
    return getSafeIntFromInput(torusToolbox.torusRingCountInput);
}

int MainWindow::getTorusToolboxRingSubdiv()
{
    return getSafeIntFromInput(torusToolbox.torusRingSubdivInput);
}

double MainWindow::getCapsuleToolboxRadius()
{
    return getSafeDoubleFromInput(capsuleToolbox.capsuleRadiusInput);
}

double MainWindow::getCapsuleToolboxCylinderHeight()
{
    return getSafeDoubleFromInput(capsuleToolbox.capsuleCylinderHeightInput);
}

int MainWindow::getCapsuleToolboxCylinderHeightSubdiv()
{
    return getSafeDoubleFromInput(capsuleToolbox.capsuleCylinderHeightSubdivInput);
}

int MainWindow::getCapsuleToolboxCylinderSubdiv()
{
    return getSafeDoubleFromInput(capsuleToolbox.capsuleCylinderSubdivInput);
}

int MainWindow::getCapsuleToolboxCapsSubdiv()
{
    return getSafeDoubleFromInput(capsuleToolbox.capsuleCapsSubdivInput);
}

void MainWindow::UpdateIcosphere()
{
    double radius = this->getIcosphereToolboxRadius();
    int subdivisions = this->getIcosphereToolboxSubdiv();

    if(radius == -1 || subdivisions == -1)
        return;//Incorrect parameters, not doing anything

    CreateIcosphereMesh(radius, subdivisions);

    UpdateGeometry();
}

void MainWindow::UpdateTorus()
{
    double innerRadius = this->getTorusToolboxInRadius();
    double outerRadius = this->getTorusToolboxOutRadius();
    int ringCount = this->getTorusToolboxRingCount();
    int ringsSubdivisions = this->getTorusToolboxRingSubdiv();

    if(innerRadius == -1 || outerRadius == -1 || ringCount == -1 || ringsSubdivisions == -1)
        return;//Incorrect parameters, not doing anything

    CreateTorusMesh(innerRadius, outerRadius, ringCount, ringsSubdivisions);

    UpdateGeometry();
}

void MainWindow::UpdateCapsule()
{
    double radius = this->getCapsuleToolboxRadius();
    double cylinderHeight = this->getCapsuleToolboxCylinderHeight();
    int cylinderHeightSubdivions = this->getCapsuleToolboxCylinderHeightSubdiv();
    int cylinderSubdivisions = this->getCapsuleToolboxCylinderSubdiv();
    int sphereHeightSubdivisions = this->getCapsuleToolboxCapsSubdiv();

    if(radius == -1 || cylinderHeight == -1 ||
       cylinderHeightSubdivions == -1 || cylinderSubdivisions == -1 ||
       sphereHeightSubdivisions == -1)
        return;//Incorrect parameters, not doing anything

    CreateCapsuleMesh(radius, cylinderHeight, cylinderHeightSubdivions, cylinderSubdivisions, sphereHeightSubdivisions);

    UpdateGeometry();
}
