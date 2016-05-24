TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis_set.cpp \
    ccdnaive.cpp \
    ccdintermediates.cpp \
    ccdblocks.cpp \
    mbptnaive.cpp \
    fci.cpp

HEADERS += \
    basis_set.h \
    ccdnaive.h \
    ccdintermediates.h \
    ccdblocks.h \
    mbptnaive.h \
    fci.h

LIBS += -llapack -lblas -larmadillo
