#ifndef __Qte__
#define __Qte__

#include <QtWidgets/qmainwindow.h>
#include "realtime.h"
#include "meshcolor.h"
#include "ui_icosphereToolbox.h"

QT_BEGIN_NAMESPACE
namespace Ui { class Assets; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
private:
    Ui::Assets* uiw;                    //!< Interface
    Ui::IcosphereToolbox icoToolbox;    //!< Icosphere toolbox widget

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
    void SetupIcosphereToolbox();
    double getIcosphereToolboxRadius();
    int getIcosphereToolboxSubdiv();

public slots:
    void editingSceneLeft(const Ray&);
    void editingSceneRight(const Ray&);
    void BoxMeshExample();
    void SphereImplicitExample();
    void DisplayIcosphere();
    void ResetCamera();
    void UpdateMaterial();
    void UpdateIcosphere();
};

#endif
