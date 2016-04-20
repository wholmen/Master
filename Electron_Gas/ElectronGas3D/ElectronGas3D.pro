TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis_set.cpp \
    basis_set2d.cpp \
    naivesolvers.cpp \
    ccdintermediate.cpp \
    ccdblocks.cpp

HEADERS += \
    basis_set.h \
    basis_set2d.h \
    naivesolvers.h \
    ccdintermediate.h \
    ccdblocks.h

LIBS += -llapack -lblas -larmadillo
