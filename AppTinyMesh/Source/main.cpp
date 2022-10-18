#include "qtemainwindow.h"
#include <QtWidgets/qapplication.h>

#include "benchmarks.h"

int main(int argc, char *argv[])
{
//    int iterations = 1000;

//    Benchmarks benchmarks;

//    benchmarks.BenchmarkIcosphere(iterations, 11);
//    benchmarks.BenchmarkTorus(iterations, 600, 20, 20, false);
//    benchmarks.BenchmarkTorus(600, 600, 20, 20, true);
//    benchmarks.BenchmarkCapsule(7500, 300);
//    benchmarks.BenchmarkCylinder(500, 750);

//    return 0;

	QApplication app(argc, argv);

	MainWindow mainWin;
	mainWin.showMaximized();

	return app.exec();
}
