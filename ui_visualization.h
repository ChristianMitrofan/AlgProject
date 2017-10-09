/********************************************************************************
** Form generated from reading UI file 'visualization.ui'
**
** Created by: Qt User Interface Compiler version 5.2.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_VISUALIZATION_H
#define UI_VISUALIZATION_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Visualization
{
public:
    QAction *actionOpen;
    QAction *actionExit;
    QAction *actionInsert;
    QWidget *centralwidget;
    QTextBrowser *equationsTxt;
    QLabel *equationsLbl;
    QTextBrowser *outputTxt;
    QLabel *outputLbl;
    QPushButton *pushButton;
    QPushButton *pushButton_2;
    QTextBrowser *rootsTxt;
    QLabel *rootsLbl;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuGnupoints;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *Visualization)
    {
        if (Visualization->objectName().isEmpty())
            Visualization->setObjectName(QStringLiteral("Visualization"));
        Visualization->resize(1024, 768);
        actionOpen = new QAction(Visualization);
        actionOpen->setObjectName(QStringLiteral("actionOpen"));
        actionOpen->setVisible(true);
        actionExit = new QAction(Visualization);
        actionExit->setObjectName(QStringLiteral("actionExit"));
        actionInsert = new QAction(Visualization);
        actionInsert->setObjectName(QStringLiteral("actionInsert"));
        centralwidget = new QWidget(Visualization);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        equationsTxt = new QTextBrowser(centralwidget);
        equationsTxt->setObjectName(QStringLiteral("equationsTxt"));
        equationsTxt->setGeometry(QRect(9, 77, 1001, 81));
        equationsLbl = new QLabel(centralwidget);
        equationsLbl->setObjectName(QStringLiteral("equationsLbl"));
        equationsLbl->setGeometry(QRect(10, 40, 71, 17));
        outputTxt = new QTextBrowser(centralwidget);
        outputTxt->setObjectName(QStringLiteral("outputTxt"));
        outputTxt->setGeometry(QRect(20, 260, 1001, 271));
        outputLbl = new QLabel(centralwidget);
        outputLbl->setObjectName(QStringLiteral("outputLbl"));
        outputLbl->setGeometry(QRect(20, 230, 71, 17));
        pushButton = new QPushButton(centralwidget);
        pushButton->setObjectName(QStringLiteral("pushButton"));
        pushButton->setGeometry(QRect(750, 20, 251, 27));
        pushButton_2 = new QPushButton(centralwidget);
        pushButton_2->setObjectName(QStringLiteral("pushButton_2"));
        pushButton_2->setGeometry(QRect(230, 170, 99, 27));
        rootsTxt = new QTextBrowser(centralwidget);
        rootsTxt->setObjectName(QStringLiteral("rootsTxt"));
        rootsTxt->setGeometry(QRect(20, 570, 1001, 141));
        rootsLbl = new QLabel(centralwidget);
        rootsLbl->setObjectName(QStringLiteral("rootsLbl"));
        rootsLbl->setGeometry(QRect(20, 540, 111, 17));
        Visualization->setCentralWidget(centralwidget);
        menubar = new QMenuBar(Visualization);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 1024, 25));
        menubar->setDefaultUp(true);
        menubar->setNativeMenuBar(false);
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        menuGnupoints = new QMenu(menubar);
        menuGnupoints->setObjectName(QStringLiteral("menuGnupoints"));
        Visualization->setMenuBar(menubar);
        statusbar = new QStatusBar(Visualization);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        Visualization->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuGnupoints->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionExit);
        menuGnupoints->addSeparator();
        menuGnupoints->addAction(actionInsert);

        retranslateUi(Visualization);

        QMetaObject::connectSlotsByName(Visualization);
    } // setupUi

    void retranslateUi(QMainWindow *Visualization)
    {
        Visualization->setWindowTitle(QApplication::translate("Visualization", "Visualization App", 0));
        actionOpen->setText(QApplication::translate("Visualization", "Open", 0));
        actionOpen->setShortcut(QApplication::translate("Visualization", "Alt+F", 0));
        actionExit->setText(QApplication::translate("Visualization", "Exit", 0));
        actionExit->setShortcut(QApplication::translate("Visualization", "Alt+X", 0));
        actionInsert->setText(QApplication::translate("Visualization", "Insert", 0));
        equationsLbl->setText(QApplication::translate("Visualization", "Equations", 0));
        outputLbl->setText(QApplication::translate("Visualization", "Output", 0));
        pushButton->setText(QApplication::translate("Visualization", "Plot Algebraic Curves", 0));
        pushButton_2->setText(QApplication::translate("Visualization", "Rm Pol/mial", 0));
        rootsLbl->setText(QApplication::translate("Visualization", " Roots (x , y)", 0));
        menuFile->setTitle(QApplication::translate("Visualization", "File", 0));
        menuGnupoints->setTitle(QApplication::translate("Visualization", "Gnupoints", 0));
    } // retranslateUi

};

namespace Ui {
    class Visualization: public Ui_Visualization {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_VISUALIZATION_H
