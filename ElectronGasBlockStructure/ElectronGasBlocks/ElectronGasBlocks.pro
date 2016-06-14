TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    solver.cpp \
    block.cpp \
    electronbasis.cpp \
    ../../Solvers/Solvers/mbptnaive.cpp \
    ../../Solvers/Solvers/ccdnaive.cpp \
    ../../Solvers/Solvers/ccdintermediates.cpp \
    ../../Pairing_Model/PairingModel/pairingbasis.cpp \
    ../../NuclearMatter/NuclearMatter/nuclearbasis.cpp

HEADERS += \
    solver.h \
    block.h \
    electronbasis.h \
    ../../Solvers/Solvers/mbptnaive.h \
    ../../Solvers/Solvers/ccdnaive.h \
    ../../Solvers/Solvers/ccdintermediates.h \
    ../../Pairing_Model/PairingModel/pairingbasis.h \
    ../../NuclearMatter/NuclearMatter/nuclearbasis.h


LIBS += -llapack -lblas -larmadillo
