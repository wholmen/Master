TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    fci.cpp \
    pairingbasis.cpp \
    ../../Solvers/Solvers/ccdintermediates.cpp \
    ../../Solvers/Solvers/ccdnaive.cpp \
    ../../Solvers/Solvers/mbptnaive.cpp \
    ../../ElectronGasBlockStructure/ElectronGasBlocks/electronbasis.cpp \
    ../../NuclearMatter/NuclearMatter/nuclearbasis.cpp

HEADERS += \
    fci.h \
    pairingbasis.h \
    ../../Solvers/Solvers/ccdintermediates.h \
    ../../Solvers/Solvers/ccdnaive.h \
    ../../Solvers/Solvers/mbptnaive.h \
    ../../ElectronGasBlockStructure/ElectronGasBlocks/electronbasis.h \
    ../../NuclearMatter/NuclearMatter/nuclearbasis.h

LIBS += -llapack -lblas -larmadillo
