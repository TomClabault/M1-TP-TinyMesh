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
    connect(uiw->cylinderButton, SIGNAL(clicked()), this, SLOT(DisplayCylinder()));
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
    boxMesh.accessibility(cols, 2, 15);

    //for (size_t i = 0; i < cols.size(); i++)
        //cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

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
    implicitMesh.accessibility(cols, 1, 15);

//    for (size_t i = 0; i < cols.size(); i++)
//        cols[i] = Color(0.8, 0.8, 0.8);

    meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
    UpdateGeometry();
}

bool merging = true;

void MainWindow::CreateIcosphereMesh(double radius, int subdivisions)
{
    Mesh icosphereMesh = Mesh(Icosphere(radius, subdivisions));

    std::vector<Color> cols;
    cols.resize(icosphereMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(icosphereMesh, cols, icosphereMesh.VertexIndexes());
}

/**
 *  Code pour la génération de l'union de primitives avec transformations / déformations utilisée dans le rapport
    Mesh torusMesh = Mesh(Torus(innerRadius, outerRadius, ringCount, ringsSubdivisions));

    Mesh icosphereMesh(Icosphere(1, 4));
    icosphereMesh.SphereWarp(Sphere(2, Vector(0, 0, 1)));
    icosphereMesh.SphereWarp(Sphere(3, Vector(0, 0, -2)));
    torusMesh.Merge(icosphereMesh);

    Mesh capsuleMesh = Mesh(Capsule(1, 2, 5, 10, 10));
    capsuleMesh.Rotate(Matrix::RotationX(45));
    capsuleMesh.Translate(Vector(0, 2, 2));
    torusMesh.Merge(capsuleMesh);

    Mesh capsuleMesh2 = Mesh(Capsule(1, 2, 5, 10, 10));
    capsuleMesh2.Rotate(Matrix::RotationX(-45));
    capsuleMesh2.Translate(Vector(0, -3, 0.25));
    torusMesh.Merge(capsuleMesh2);

    Mesh cylinderMesh = Mesh(Cylinder(1, 2, 4, 15));
    cylinderMesh.Translate(Vector(-3, -1, 1));
    torusMesh.Merge(cylinderMesh);
*/

void MainWindow::CreateTorusMesh(double innerRadius, double outerRadius, int ringCount, int ringsSubdivisions)
{
    Mesh torusMesh = Mesh(Torus(innerRadius, outerRadius, ringCount, ringsSubdivisions));

    std::vector<Color> cols;
    cols.resize(torusMesh.Vertexes());

    std::srand(2);
    torusMesh.accessibility(cols, 1, 10, 1);

    //TODO remove
    /*bool found = false;
    int seed = 0;

    while (!found)
    {
        std::srand(seed);
        torusMesh.accessibility(cols, 1, 1, 50);
        for (Color color : cols)
        {
            if ((color[0] != 1.0 || color[1] != 1.0 || color[2] != 1.0) && !(color[0] == 1.0 && color[1] == 0.0 && color[2] == 0.0))
            {
                found = true;

                std::cout << "Seed: " << seed << std::endl;

                std::cout << "Colors: [" << color[0] << ", " << color[1] << ", " << color[2] << "]" << std::endl;

                break;
            }
        }

        seed++;
    }*/
    

    meshColor = MeshColor(torusMesh, cols, torusMesh.VertexIndexes());
}

void MainWindow::CreateCapsuleMesh(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions)
{
    Mesh capsuleMesh = Mesh(Capsule(radius, cylinderHeight, cylinderHeightSubdivions, cylinderSubdivisions, sphereHeightSubdivisions));

    std::vector<Color> cols;
    cols.resize(capsuleMesh.Vertexes());
    capsuleMesh.accessibility(cols, 1, 15);
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(capsuleMesh, cols, capsuleMesh.VertexIndexes());
}

void MainWindow::CreateCylinderMesh(double radius, double height, int heightSubdivisions, int cylinderSubdivisions)
{
    Mesh humanMesh;
    humanMesh.Load("Skull.obj");

    std::vector<Color> cols;
    cols.resize(humanMesh.Vertexes());

    humanMesh.accessibility(cols, 1, 10);

    //Mesh cylinderMesh = Mesh(Cylinder(radius, height, heightSubdivisions, cylinderSubdivisions));

    //std::vector<Color> cols;
    //cols.resize(cylinderMesh.Vertexes());
    //for (size_t i = 0; i < cols.size(); i++)
//        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    //meshColor = MeshColor(cylinderMesh, cols, cylinderMesh.VertexIndexes());
    meshColor = MeshColor(humanMesh, cols, humanMesh.VertexIndexes());
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

void MainWindow::SetupCylinderToolbox()
{
    uiw->toolboxGroupBox->setVisible(true);

    delete toolboxWidget;//Deleting the previous widget
    toolboxWidget = new QWidget;
    cylinderToolbox.setupUi(toolboxWidget);

    QVBoxLayout vboxLayout(uiw->toolboxGroupBox);
    vboxLayout.addWidget(toolboxWidget);

    //Default settings for the icosphere at its creation
    cylinderToolbox.cylinderRadiusInput->setText("1.0");
    cylinderToolbox.cylinderHeightInput->setText("2.0");
    cylinderToolbox.cylinderHeightSubdivInput->setText("5");
    cylinderToolbox.cylinderSubdivInput->setText("10");

    connect(cylinderToolbox.cylinderRadiusInput, SIGNAL(returnPressed()), this, SLOT(UpdateCylinder()));
    connect(cylinderToolbox.cylinderHeightInput, SIGNAL(returnPressed()), this, SLOT(UpdateCylinder()));
    connect(cylinderToolbox.cylinderHeightSubdivInput, SIGNAL(returnPressed()), this, SLOT(UpdateCylinder()));
    connect(cylinderToolbox.cylinderSubdivInput, SIGNAL(returnPressed()), this, SLOT(UpdateCylinder()));
    connect(cylinderToolbox.applyCylinderToolboxButton, SIGNAL(clicked()), this, SLOT(UpdateCylinder()));
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

void MainWindow::DisplayCylinder()
{
    CreateCylinderMesh(1.0, 2.0, 5, 10);

    UpdateGeometry();

    SetupCylinderToolbox();
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

void MainWindow::UpdateIcosphere()
{
    double radius = getSafeDoubleFromInput(icosphereToolbox.icosphereRadiusInput);
    int subdivisions = getSafeIntFromInput(icosphereToolbox.icosphereSubdivisionInput);

    if(radius == -1 || subdivisions == -1)
        return;//Incorrect parameters, not doing anything

    CreateIcosphereMesh(radius, subdivisions);

    UpdateGeometry();
}

void MainWindow::UpdateTorus()
{
    double innerRadius = getSafeDoubleFromInput(torusToolbox.torusInRadiusInput);
    double outerRadius = getSafeDoubleFromInput(torusToolbox.torusOutRadiusInput);
    int ringCount = getSafeIntFromInput(torusToolbox.torusRingCountInput);
    int ringsSubdivisions = getSafeIntFromInput(torusToolbox.torusRingSubdivInput);

    if(innerRadius == -1 || outerRadius == -1 || ringCount == -1 || ringsSubdivisions == -1)
        return;//Incorrect parameters, not doing anything

    CreateTorusMesh(innerRadius, outerRadius, ringCount, ringsSubdivisions);

    UpdateGeometry();
}

void MainWindow::UpdateCapsule()
{
    double radius = getSafeDoubleFromInput(capsuleToolbox.capsuleRadiusInput);
    double cylinderHeight = getSafeDoubleFromInput(capsuleToolbox.capsuleCylinderHeightInput);
    int cylinderHeightSubdivions = getSafeIntFromInput(capsuleToolbox.capsuleCylinderHeightSubdivInput);
    int cylinderSubdivisions = getSafeIntFromInput(capsuleToolbox.capsuleCylinderSubdivInput);
    int capsSubdivisions = getSafeIntFromInput(capsuleToolbox.capsuleCapsSubdivInput);

    if(radius == -1 || cylinderHeight == -1 ||
            cylinderHeightSubdivions == -1 || cylinderSubdivisions == -1 ||
            capsSubdivisions == -1)
        return;//Incorrect parameters, not doing anything

    CreateCapsuleMesh(radius, cylinderHeight, cylinderHeightSubdivions, cylinderSubdivisions, capsSubdivisions);

    UpdateGeometry();
}

void MainWindow::UpdateCylinder()
{
    double radius = getSafeDoubleFromInput(cylinderToolbox.cylinderRadiusInput);
    double height = getSafeDoubleFromInput(cylinderToolbox.cylinderHeightInput);
    int heightSubdivions = getSafeDoubleFromInput(cylinderToolbox.cylinderHeightSubdivInput);
    int subdivisions = getSafeDoubleFromInput(cylinderToolbox.cylinderSubdivInput);

    if(radius == -1 || height == -1 ||
            heightSubdivions == -1 || subdivisions == -1)
        return;//Incorrect parameters, not doing anything

    CreateCylinderMesh(radius, height, heightSubdivions, subdivisions);

    UpdateGeometry();
}
