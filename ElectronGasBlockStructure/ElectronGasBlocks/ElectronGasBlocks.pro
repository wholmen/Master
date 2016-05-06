TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis_set.cpp \
    solver.cpp \
    block.cpp

HEADERS += \
    basis_set.h \
    solver.h \
    block.h


LIBS += -llapack -lblas -larmadillo
