TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis_set.cpp \
    naivesolvers.cpp

HEADERS += \
    basis_set.h \
    naivesolvers.h

LIBS += -llapack -lblas -larmadillo
