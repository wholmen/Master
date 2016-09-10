TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    solver.cpp \
    electronbasis.cpp \
    block.cpp \
    ../../Pairing_Model/PairingModel/pairingbasis.cpp \
    ../../Solvers/Solvers/ccdintermediates.cpp \
    ../../Solvers/Solvers/ccdnaive.cpp \
    ../../NuclearMatter/NuclearMatter/nuclearbasis.cpp

HEADERS += \
    solver.h \
    electronbasis.h \
    block.h \
    ../../NuclearMatter/NuclearMatter/nuclearbasis.h \
    ../../Pairing_Model/PairingModel/pairingbasis.h \
    ../../Solvers/Solvers/ccdintermediates.h \
    ../../Solvers/Solvers/ccdnaive.h

LIBS += -llapack -lblas -larmadillo


QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
LIBS += -fopenmp
