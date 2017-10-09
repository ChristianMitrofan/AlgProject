//#include <QApplication>
#include "visualization.h"
#include "Polynomials.h"

int main(int argc, char* argv[])
{
    setlocale(LC_NUMERIC,"en_US");
	QApplication app(argc, argv);
    QCoreApplication::setAttribute(Qt::AA_DontUseNativeMenuBar);
    Visualization v;
    v.show();

    return app.exec();
}
