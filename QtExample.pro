QT = core gui widgets network printsupport svg

TARGET = qtgnuplotlib-example
TEMPLATE = app
LIBS += -lQtGnuplot \
            /usr/lib/libblas.so \
            /usr/lib/liblapack.so \
            /usr/lib/liblapacke.so \

INCLUDEPATH += ./eigen-eigen-b30b87236a1b


CONFIG += debug

SOURCES = qtgnuplotlib-example.cpp \
    visualization.cpp \
    Calculations.cpp \
    Companion.cpp \
    GenericMatrix.cpp \
    Lmatrices.cpp \
    Multiply.cpp \
    points.cpp \
    Polynomials.cpp \
    Solution.cpp \
    Sylvester.cpp \
    UnivariateSolving.cpp \
    VectorCompanion.cpp

FORMS += \
    visualization.ui

HEADERS += \
    visualization.h \
    Calculations.h \
    Companion.h \
    GenericMatrix.h \
    Lmatrices.h \
    Multiply.h \
    points.h \
    Polynomials.h \
    Solution.h \
    Sylvester.h \
    UnivariateSolving.h \
    VectorCompanion.h

