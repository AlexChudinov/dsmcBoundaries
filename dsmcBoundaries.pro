#-------------------------------------------------
#
# Project created by QtCreator 2020-09-16T15:21:52
#
#-------------------------------------------------

QT       -= core gui

TARGET = dsmcFoam
TEMPLATE = app

DEFINES += DSMCPROJECT_LIBRARY

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    BinaryCollisionModel.C \
    BinaryCollisionModelNew.C \
    defineDSMCParcel.C \
    DSMCCloud.C \
    DSMCCloudName.C \
    dsmcFoam.C \
    dsmcParcel.C \
    DSMCParcel.C \
    DSMCParcelIO.C \
    FreeStream.C \
    InflowBoundaryModel.C \
    InflowBoundaryModelNew.C \
    LarsenBorgnakkeVariableHardSphere.C \
    makeDSMCParcelBinaryCollisionModels.C \
    makeDSMCParcelInflowBoundaryModels.C \
    makeDSMCParcelWallInteractionModels.C \
    MaxwellianThermal.C \
    MixedDiffuseSpecular.C \
    NoBinaryCollision.C \
    NoInflow.C \
    SpecularReflection.C \
    VariableHardSphere.C \
    WallInteractionModel.C \
    WallInteractionModelNew.C

HEADERS += \
    BinaryCollisionModel.H \
    dsmcCloud.H \
    DSMCCloud.H \
    dsmcParcel.H \
    DSMCParcel.H \
    DSMCParcelI.H \
    FreeStream.H \
    InflowBoundaryModel.H \
    LarsenBorgnakkeVariableHardSphere.H \
    MaxwellianThermal.H \
    MixedDiffuseSpecular.H \
    NoBinaryCollision.H \
    NoInflow.H \
    SpecularReflection.H \
    VariableHardSphere.H \
    WallInteractionModel.H \
    inputstream.h \
    createFields.H

unix {
    target.path = /usr/lib
    INSTALLS += target
}

DEFINES += \
    linux64 \
    WM_ARCH_OPTION=64 \
    WM_DP \
    WM_LABEL_SIZE=32 \
    NoRepository \

INCLUDEPATH += \
    . \
    /opt/openfoam-dev/src/finiteVolume/lnInclude \
    /opt/openfoam-dev/src/lagrangian/basic/lnInclude \
    /opt/openfoam-dev/src/meshTools/lnInclude \
    /opt/openfoam-dev/src/OpenFOAM/lnInclude \
    /opt/openfoam-dev/src/OSspecific/POSIX/lnInclude

QMAKE_CXXFLAGS += \
    -std=c++14 \
    -m64 \
    -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter \
    -Wno-invalid-offsetof -Wno-attributes -O3  -ftemplate-depth-100 \
    -fPIC -c

QMAKE_LIBDIR += /opt/openfoam-dev/platforms/linux64GccDPInt32Opt/lib \
    /opt/openfoam-dev/platforms/linux64GccDPInt32Opt/lib/dummy

LIBS += -lmeshTools -lfiniteVolume -llagrangian -lDSMC \
    -lOpenFOAM -ldl -ltriSurface -lsurfMesh -lPstream -lfileFormats
