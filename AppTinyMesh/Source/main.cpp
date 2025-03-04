#include "analyticApproximations.h"
#include "benchmarks.h"
#include "BVH.h"

#include "qtemainwindow.h"
#include <QtWidgets/qapplication.h>

int main(int argc, char *argv[])
{
//    AnalyticSphere::tests();
//    AnalyticCylinder::tests();
//    BVH::BVHTests();
//    return 0;

    //Benchmarks::BenchmarkBVHIntersectionCount();
    //return 0;
#ifdef VISUAL_STUDIO_EDITOR
	//Opening a console so that we can get cout output
	AllocConsole();
	freopen("CONOUT$", "w", stdout);
	freopen("CONOUT$", "w", stderr);
#endif

	QApplication app(argc, argv);

	MainWindow mainWin;
	mainWin.showMaximized();

	return app.exec();
}
