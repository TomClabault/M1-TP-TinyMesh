QT += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
greaterThan(QT_MAJOR_VERSION, 5): QT += openglwidgets

CONFIG += c++11

INCLUDEPATH += AppTinyMesh/Include
INCLUDEPATH += $$(GLEW_DIR)
INCLUDEPATH += $$(OUT_PWD)

VPATH += AppTinyMesh

SOURCES += \
    AppTinyMesh/Source/analyticApproximations.cpp \
    AppTinyMesh/Source/benchmarks.cpp \
    AppTinyMesh/Source/box.cpp \
    AppTinyMesh/Source/color.cpp \
    AppTinyMesh/Source/evector.cpp \
    AppTinyMesh/Source/implicits.cpp \
    AppTinyMesh/Source/main.cpp \
    AppTinyMesh/Source/camera.cpp \
    AppTinyMesh/Source/mathematics.cpp \
    AppTinyMesh/Source/matrix.cpp \
    AppTinyMesh/Source/mesh.cpp \
    AppTinyMesh/Source/meshcolor.cpp \
    AppTinyMesh/Source/mesh-widget.cpp \
    AppTinyMesh/Source/qtemainwindow.cpp \
    AppTinyMesh/Source/ray.cpp \
    AppTinyMesh/Source/simpleMeshes.cpp \
    AppTinyMesh/Source/shader-api.cpp \
    AppTinyMesh/Source/sphere.cpp \
    AppTinyMesh/Source/triangle.cpp \

HEADERS += \
    AppTinyMesh/Include/analyticApproximations.h \
    AppTinyMesh/Include/benchmarks.h \
    AppTinyMesh/Include/box.h \
    AppTinyMesh/Include/camera.h \
    AppTinyMesh/Include/color.h \
    AppTinyMesh/Include/implicits.h \
    AppTinyMesh/Include/mathematics.h \
    AppTinyMesh/Include/matrix.h \
    AppTinyMesh/Include/mesh.h \
    AppTinyMesh/Include/meshcolor.h \
    AppTinyMesh/Include/qtemainwindow.h \
    AppTinyMesh/Include/realtime.h \
    AppTinyMesh/Include/simpleMeshes.h \
    AppTinyMesh/Include/shader-api.h \
    AppTinyMesh/Include/sphere.h \

FORMS += \
    AppTinyMesh/UI/icosphereToolbox.ui \
    AppTinyMesh/UI/torusToolbox.ui \
    AppTinyMesh/UI/capsuleToolbox.ui \
    AppTinyMesh/UI/cylinderToolbox.ui \
    AppTinyMesh/UI/interface.ui

win32 {
    LIBS += -L$$(GLEW_DIR) -lglew32
    LIBS += -lopengl32 -lglu32
}
unix:!macx {
    LIBS += -lGLEW -lGL -lGLU
}
macx {
    LIBS += -lGLEW -lGL -lGLU
}

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
