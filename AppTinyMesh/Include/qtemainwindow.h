#ifndef __Qte__
#define __Qte__

#include <QtWidgets/qmainwindow.h>
#include "realtime.h"
#include "meshcolor.h"
#include "ui_icosphereToolbox.h"
#include "ui_torusToolbox.h"
#include "ui_capsuleToolbox.h"
#include "ui_cylinderToolbox.h"

QT_BEGIN_NAMESPACE
namespace Ui { class Assets; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
private:
    Ui::Assets* uiw;                    //!< Interface
    Ui::IcosphereToolbox icosphereToolbox;    //!< Icosphere toolbox widget
    Ui::TorusToolbox torusToolbox;    //!< Torus toolbox widget
    Ui::CapsuleToolbox capsuleToolbox;    //!< Capsule toolbox widget
    Ui::CylinderToolbox cylinderToolbox;    //!< Capsule toolbox widget

    QWidget* toolboxWidget = nullptr;   //!< Widget that allows the modification of
    //!< the object in the current meshWidget
    //! (subdivisions, radius, ...)
    MeshWidget* meshWidget;   //!< Viewer
    MeshColor meshColor;		//!< Mesh.

public:
    MainWindow();
    ~MainWindow();
    void CreateActions();
    void UpdateGeometry();

    void CreateIcosphereMesh(double radius, int subdivisions);
    void CreateTorusMesh(double innerRadius, double outerRadius, int ringCount, int ringsSubdivisions);
    void CreateCapsuleMesh(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions);
    void CreateCylinderMesh(double radius, double height, int heightSubdivisions, int cylinderSubdivisions);

    void SetupIcosphereToolbox();
    void SetupTorusToolbox();
    void SetupCapsuleToolbox();
    void SetupCylinderToolbox();

public slots:
    void editingSceneLeft(const Ray&);
    void editingSceneRight(const Ray&);
    void BoxMeshExample();
    void SphereImplicitExample();

    void DisplayIcosphere();
    void DisplayTorus();
    void DisplayCapsule();
    void DisplayCylinder();

    void ResetCamera();
    void UpdateMaterial();

    void UpdateIcosphere();
    void UpdateTorus();
    void UpdateCapsule();
    void UpdateCylinder();
};

#endif
