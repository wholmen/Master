TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis_set.cpp \
    ccdintermediates.cpp \
    ccdblocks.cpp

HEADERS += \
    basis_set.h \
    ccdintermediates.h \
    ccdblocks.h

