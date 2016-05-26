TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis_set.cpp \
    solver.cpp \
    block.cpp \
    ccdintermediates.cpp

HEADERS += \
    basis_set.h \
    solver.h \
    block.h \
    ccdintermediates.h


LIBS += -llapack -lblas -larmadillo
