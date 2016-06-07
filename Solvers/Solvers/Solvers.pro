TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    ccdnaive.cpp \
    ccdintermediates.cpp \
    mbptnaive.cpp \
    ../../Pairing_Model/PairingModel/pairingbasis.cpp \
    ../../ElectronGasBlockStructure/ElectronGasBlocks/electronbasis.cpp

HEADERS += \
    ccdnaive.h \
    ccdintermediates.h \
    mbptnaive.h \
    ../../Pairing_Model/PairingModel/pairingbasis.h \
    ../../ElectronGasBlockStructure/ElectronGasBlocks/electronbasis.h


LIBS += -llapack -lblas -larmadillo
