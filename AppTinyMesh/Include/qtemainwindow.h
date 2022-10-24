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

    /*!
     * \brief Gets and stores the ambient occlusion parameters
     * of the input of the UI in the given parameters.
     * If any invalid value is found in the inputs of the UI,
     * all the parameters will be set to default values
     *
     * \param AORadius [out] The radius of the ambient
     * occlusion will be stored in this parameter.
     * Will be set to 1 if an invalid value if found in
     * the inputs of the UI
     * \param AOSamples [out] The numbers of samples used
     * to compute the ambient occlusion will be stored in
     * this parameter. Will be set to 10 if an invalid
     * value if found in the inputs of the UI
     * \param AOStrength [out] The strength of the ambient
     * occlusion will be stored in this parameter.
     * Will be set to 1 if an invalid value if found in
     * the inputs of the UI
     */
    void GetAOParameters(double& AORadius, int& AOSamples, double& AOStrength);

    /*!
     * \brief Checks the state of the ambient occlusion
     * checkbox of the UI and computes the colors of the
     * vertices of the mesh accordingly. If the checkbox
     * is checked, the ambient occlusion of the mesh will
     * be computed using the parameters of the ambient
     * occlusion found in the inputs of the UI.
     * If the ambient occlusion checkbox isn't checked,
     * the color of all the vertices of the mesh will be
     * set to white.
     *
     * \param mesh The mesh whose ambient occlusion will
     * be computed if the ambient occlusion checkbox is checked.
     * \param meshColors The resulting color after the
     * computation of the ambient occlusion
     * (or the set operation to all white)
     */
    void HandleAO(Mesh mesh, std::vector<Color>& meshColors);

    void CreateIcosphereMesh(double radius, int subdivisions);
    void CreateTorusMesh(double innerRadius, double outerRadius, int ringCount, int ringsSubdivisions);
    void CreateCapsuleMesh(double radius, double cylinderHeight, int cylinderHeightSubdivions, int cylinderSubdivisions, int sphereHeightSubdivisions);
    void CreateCylinderMesh(double radius, double height, int heightSubdivisions, int cylinderSubdivisions);

    void LoadObjMesh(QString objFilePath, double occlusionRadius, int occlusionSamples, double occlusionStrength);

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

    void DisplayObjMesh();

    void ResetCamera();
    void UpdateMaterial();
    void UpdateAO();

    void UpdateIcosphere();
    void UpdateTorus();
    void UpdateCapsule();
    void UpdateCylinder();
};

#endif
