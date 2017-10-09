#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <QMainWindow>
#include <QtGnuplot/QtGnuplotWidget.h>
#include <QtGnuplot/QtGnuplotInstance.h>
#include <GenericMatrix.h>
//#include <GenericMatrix.cpp>
#include <Eigen/Core>


namespace Ui {
class Visualization;
}

class Visualization : public QMainWindow
{
    Q_OBJECT

public:
    explicit Visualization(QWidget *parent = 0);
    ~Visualization();

private:
    QtGnuplotWidget *widget;
    QtGnuplotInstance *instance;
    QPushButton *solve;
    int points;
    int d;
    double ** point;
    bool solved;
    void print_results(GenericMatrix<double *> *syl, SylvPol<GenericMatrix<double> > *sylp, double k, char hidden);

protected:
    bool eventFilter(QObject *obj, QEvent *event);

private slots:

    void replot_points();

    void solve_clicked();

    void on_actionExit_triggered();

    void on_actionOpen_triggered();

    void on_pushButton_clicked();

    void on_actionInsert_triggered();

    void on_equationsTxt_textChanged();

    void pointsClosed();

    void on_pushButton_2_clicked();

private:
    Ui::Visualization *ui;
};

#endif // VISUALIZATION_H
