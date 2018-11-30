MOC_DIR = ./moc
OBJECTS_DIR = ./obj

HEADERS       = window.h \
                mainwindow.h \
    BasicPoint.h \
    Voxel.h \
    VoxelGrid.h \
    TetMeshCreator.h \
    Tetrahedron.h
SOURCES       = main.cpp \
                window.cpp \
                mainwindow.cpp \
    VoxelGrid.cpp \
    TetMeshCreator.cpp

QT           += widgets
