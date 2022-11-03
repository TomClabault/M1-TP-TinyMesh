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
    connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));

    connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw->icosphereButton, SIGNAL(clicked()), this, SLOT(DisplayIcosphere()));
    connect(uiw->torusButton, SIGNAL(clicked()), this, SLOT(DisplayTorus()));
    connect(uiw->capsuleButton, SIGNAL(clicked()), this, SLOT(DisplayCapsule()));
    connect(uiw->cylinderButton, SIGNAL(clicked()), this, SLOT(DisplayCylinder()));

    connect(uiw->lotusFlowerButton, SIGNAL(clicked()), this, SLOT(DisplayObjMesh()));
    connect(uiw->lowpolyTreeButton, SIGNAL(clicked()), this, SLOT(DisplayObjMesh()));

    connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

    connect(uiw->AOCheckbox, SIGNAL(clicked()), this, SLOT(UpdateAO()));
    connect(uiw->AORadiusInput, SIGNAL(returnPressed()), this, SLOT(UpdateAO()));
    connect(uiw->AOSamplesInput, SIGNAL(returnPressed()), this, SLOT(UpdateAO()));
    connect(uiw->AOStrengthInput, SIGNAL(returnPressed()), this, SLOT(UpdateAO()));

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

double getSafeDoubleFromInput(const QLineEdit* input)
{
    bool ok;

    double value = input->text().toDouble(&ok);

    if(!ok)
        return -1;
    else
        return value;
}

int getSafeIntFromInput(const QLineEdit* input)
{
    bool ok;

    int value = input->text().toInt(&ok);

    if(!ok)
        return -1;
    else
        return value;
}

void MainWindow::GetAOParameters(double& AORadius, int& AOSamples, double& AOStrength)
{
    AORadius = getSafeDoubleFromInput(uiw->AORadiusInput);
    AOSamples = getSafeIntFromInput(uiw->AOSamplesInput);
    AOStrength = getSafeDoubleFromInput(uiw->AOStrengthInput);

    if(AORadius == -1 || AOSamples == -1 || AOStrength == -1)//There were invalid values in the inputs
    {
        //Going with default values
        AORadius = 1;
        AOSamples = 10;
        AOSamples = 1;
    }
}

void MainWindow::HandleAO(Mesh& mesh, std::vector<Color>& cols)
{
    double AORadius;
    int AOSamples;
    double AOStrength;
    GetAOParameters(AORadius, AOSamples, AOStrength);

    if(uiw->AOCheckbox->isChecked())
        mesh.accessibility(cols, AORadius, AOSamples, AOStrength);
    else
    {
        for(Color& color : cols)
            color = Color(1.0);
    }
}

void MainWindow::BoxMeshExample()
{
    uiw->toolboxGroupBox->setVisible(false);

    Mesh boxMesh = Mesh(Box(1.0));

    std::vector<Color> cols;
    cols.resize(boxMesh.Vertexes());
    HandleAO(boxMesh, cols);

    meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());

    std::cout << "indexes: ";
    for(int color : meshColor.ColorIndexes())
        std::cout << color << " ";
    std::cout << std::endl;

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
    HandleAO(implicitMesh, cols);

    meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::CreateIcosphereMesh(double radius, int subdivisions)
{
//normal: Vector(0,0,1)
//ray origin: Vector(3.07167,-0.778475,-1.27495)
//ray direction: Vector(-0.521938,0.55617,0.646727)

      //Mesh icosphereMesh = Mesh(Icosphere(Vector(radius * 2, 0, 0), radius, subdivisions));
//      icosphereMesh.Merge(Mesh(Box(Vector(3.07167,-0.778475,-1.27495), 0.025)));
//      icosphereMesh.Merge(Mesh(Box(Vector(3.07167,-0.778475,-1.27495) + Vector(-0.521938,0.55617,0.646727) * 0.2, 0.025)));
    Mesh icosphereMesh = Mesh(Icosphere(radius, subdivisions));
    icosphereMesh.Merge(Mesh(Icosphere(Vector(radius * 2, 0, 0), radius, subdivisions)));
    icosphereMesh.Merge(Mesh(Icosphere(Vector(radius, 0, std::sqrt(3) * radius), radius, subdivisions)));
    //icosphereMesh.Merge(Mesh(Icosphere(Vector(radius, 2*std::sqrt(6)/3, std::sqrt(3) / 3 * radius * radius), radius, subdivisions)));
    icosphereMesh.Merge(Mesh(Icosphere(Vector(radius, (2 * std::sqrt(6) / 3) * radius, std::sqrt(3) / 3 * radius), radius, subdivisions)));

//    icosphereMesh.Merge(Mesh(Box(Vector(0.850651,-0.525731,0.0001), 0.025)));
//    icosphereMesh.Merge(Mesh(Box(Vector(0.850651,-0.525731,0.0001) + Vector(1/2.0,0, 0), 0.025)));
//    ray origin: Vector(0.850651,-0.525731,0.0001)
//    ray direction: Vector(0.268032,0.365336,0.891453)

    //TODO debug inter volumique sphere
    //TODO debug la lenteur de la BVH
    //TODO inter volumique tore
    //TODO comparaison perf inter volumique sphere et non volumique + partie dans le rapport

    //TODO remove
//    for(int i = 0; i < 5; i++)
//    {
//        Mesh boxMesh = Mesh(Box(Vector(3.07167,-0.778475,-1.27495) + (i / 5.0) * Vector(0,0,1) * 4, 0.025));
//        icosphereMesh.Merge(boxMesh);
//    }

    std::vector<Color> cols;
    cols.resize(icosphereMesh.Vertexes());
    HandleAO(icosphereMesh, cols);

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
    HandleAO(torusMesh, cols);

    meshColor = MeshColor(torusMesh, cols, torusMesh.VertexIndexes());
}

void MainWindow::CreateCapsuleMesh(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions)
{
    Mesh capsuleMesh = Mesh(Capsule(radius, cylinderHeight, cylinderHeightSubdivions, cylinderSubdivisions, sphereHeightSubdivisions));

    std::vector<Color> cols;
    cols.resize(capsuleMesh.Vertexes());
    HandleAO(capsuleMesh, cols);

    meshColor = MeshColor(capsuleMesh, cols, capsuleMesh.VertexIndexes());
}

void MainWindow::CreateCylinderMesh(double radius, double height, int heightSubdivisions, int cylinderSubdivisions)
{
    Mesh cylinderMesh = Mesh(Cylinder(Vector(0, 0, 0), radius, height, heightSubdivisions, cylinderSubdivisions));

    std::vector<Color> cols;
    cols.resize(cylinderMesh.Vertexes());
    HandleAO(cylinderMesh, cols);

    meshColor = MeshColor(cylinderMesh, cols, cylinderMesh.VertexIndexes());
}

void MainWindow::LoadObjMesh(QString objFilePath)
{
    Mesh objMesh;
    objMesh.Load(objFilePath);

    std::vector<Color> cols;
    cols.resize(objMesh.Vertexes());
    for(Color& color : cols)
        color = Color(1.0, 0.0, 0.0);

    HandleAO(objMesh, cols);

    meshColor = MeshColor(objMesh, cols, objMesh.VertexIndexes());
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

    //Default settings for the torus at its creation
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

    //Default settings for the capsule at its creation
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

    //Default settings for the cylinder at its creation
    cylinderToolbox.cylinderRadiusInput->setText("1.0");
    cylinderToolbox.cylinderHeightInput->setText("2.0");
    cylinderToolbox.cylinderHeightSubdivInput->setText("2");
    cylinderToolbox.cylinderSubdivInput->setText("4");

    connect(cylinderToolbox.cylinderRadiusInput, SIGNAL(returnPressed()), this, SLOT(UpdateCylinder()));
    connect(cylinderToolbox.cylinderHeightInput, SIGNAL(returnPressed()), this, SLOT(UpdateCylinder()));
    connect(cylinderToolbox.cylinderHeightSubdivInput, SIGNAL(returnPressed()), this, SLOT(UpdateCylinder()));
    connect(cylinderToolbox.cylinderSubdivInput, SIGNAL(returnPressed()), this, SLOT(UpdateCylinder()));
    connect(cylinderToolbox.applyCylinderToolboxButton, SIGNAL(clicked()), this, SLOT(UpdateCylinder()));
}

void MainWindow::DisplayIcosphere()
{
    SetupIcosphereToolbox();

    UpdateIcosphere();
}

void MainWindow::DisplayTorus()
{
    SetupTorusToolbox();

    UpdateTorus();
}

void MainWindow::DisplayCapsule()
{
    SetupCapsuleToolbox();

    UpdateCapsule();
}

void MainWindow::DisplayCylinder()
{
    SetupCylinderToolbox();

    UpdateCylinder();
}

void MainWindow::DisplayObjMesh()
{
    std::string signalSender = this->sender()->objectName().toStdString();

    if(signalSender == "lotusFlowerButton")
        LoadObjMesh("LotusFlowerDecimate.obj");
    else if(signalSender == "lowpolyTreeButton")
        LoadObjMesh("Lowpoly_tree.obj");

    UpdateGeometry();
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

void MainWindow::UpdateAO()
{
    if(!uiw->AOCheckbox->isChecked())
    {
        std::vector<Color>* cols = meshColor.GetColorsVector();

        //Repainting the mesh in white
        for(Color& color : *cols)
            color = Color(1.0);
    }
    else
    {
        double AORadius;
        int AOSamples;
        double AOStrength;

        GetAOParameters(AORadius, AOSamples, AOStrength);

        if(AORadius != -1 && AOSamples != -1 && AOStrength != -1)
        {
//            Mesh icospherMesh = Mesh(Icosphere(1, 1));
//            std::vector<Color> cols;
//            cols.resize(icospherMesh.Vertexes());
//            for(Color& color : cols)
//                color = Color(1.0, 0.0, 0.0);

//            meshColor = MeshColor(icospherMesh, cols, icospherMesh.VertexIndexes());

            meshColor.computeAccessibility(AORadius, AOSamples, AOStrength);

//            //TODO remove
//            MeshColor meshColor2 = MeshColor(meshColor, meshColor.GetColors(), meshColor.VertexIndexes());//TODO remove
//            meshColor = meshColor2;//TODO remove
        }
    }

    UpdateGeometry();
}

void MainWindow::ResetCamera()
{
    meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
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
