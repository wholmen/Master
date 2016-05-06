TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis_set.cpp \
    naivesolvers.cpp \
    ccdblocks2.cpp

HEADERS += \
    basis_set.h \
    naivesolvers.h \
    ccdblocks2.h

LIBS += -llapack -lblas -larmadillo
