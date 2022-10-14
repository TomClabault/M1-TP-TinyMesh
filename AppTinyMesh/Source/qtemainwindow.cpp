#include "qtemainwindow.h"
#include "simpleMeshes.h"
#include "implicits.h"
#include "ui_interface.h"
#include "ui_icosphereToolbox.h"

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

void MainWindow::SetupIcosphereToolbox()
{
    uiw->toolboxGroupBox->setVisible(true);

    delete toolboxWidget;//Deleting the previous widget
    toolboxWidget = new QWidget;
    icoToolbox.setupUi(toolboxWidget);

    QVBoxLayout vboxLayout(uiw->toolboxGroupBox);
    vboxLayout.addWidget(toolboxWidget);

    //Default settings for the icosphere at its creation
    icoToolbox.icosphereRadiusInput->setText("1.0");
    icoToolbox.icosphereSubdivisionInput->setText("1");

    connect(icoToolbox.icosphereRadiusInput, SIGNAL(returnPressed()), this, SLOT(UpdateIcosphere()));
    connect(icoToolbox.icosphereSubdivisionInput, SIGNAL(returnPressed()), this, SLOT(UpdateIcosphere()));
    connect(icoToolbox.applyIcosphereToolboxButton, SIGNAL(clicked()), this, SLOT(UpdateIcosphere()));
}

void MainWindow::DisplayIcosphere()
{
    CreateIcosphereMesh(1.0, 1);

    UpdateGeometry();

    SetupIcosphereToolbox();
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

double MainWindow::getIcosphereToolboxRadius()
{
    double radius;

    QString icoToolboxRadiusText = icoToolbox.icosphereRadiusInput->text();
    try
    {
        radius = std::stod(icoToolboxRadiusText.toStdString());
    }
    catch (std::invalid_argument e)
    {
        return -1;
    }
    catch (std::out_of_range e)
    {
        return -1;
    }

    return radius;
}

int MainWindow::getIcosphereToolboxSubdiv()
{
    int subdivisions;

    QString icoToolboxSubdivText = icoToolbox.icosphereSubdivisionInput->text();
    try
    {
        subdivisions = std::stoi(icoToolboxSubdivText.toStdString());
    }
    catch (std::invalid_argument e)
    {
        return -1;
    }
    catch (std::out_of_range e)
    {
        return -1;
    }

    return subdivisions;
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
