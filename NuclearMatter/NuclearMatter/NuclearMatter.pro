TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ../../Solvers/Solvers/ccdintermediates.cpp \
    ../../Solvers/Solvers/ccdnaive.cpp \
    block.cpp \
    nuclearbasis.cpp \
    ../../Pairing_Model/PairingModel/pairingbasis.cpp \
    ../../ElectronGasBlockStructure/ElectronGasBlocks/electronbasis.cpp \
    solver.cpp

HEADERS += \
    ../../Solvers/Solvers/ccdintermediates.h \
    ../../Solvers/Solvers/ccdnaive.h \
    block.h \
    nuclearbasis.h \
    ../../Pairing_Model/PairingModel/pairingbasis.h \
    ../../ElectronGasBlockStructure/ElectronGasBlocks/electronbasis.h \
    solver.h

LIBS += -llapack -lblas -larmadillo


QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
LIBS += -fopenmp
