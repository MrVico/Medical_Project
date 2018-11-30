#ifndef TETMESHCREATOR_H
#define TETMESHCREATOR_H

#include "BasicPoint.h"
#include "VoxelGrid.h"
#include "Tetrahedron.h"
#include <QStatusBar>

class TetMeshCreator
{
public:
    TetMeshCreator():dataLoaded(false){
    }

    const BasicPoint & getBBMin(){ return BBMin; }
    const BasicPoint & getBBMax(){ return BBMax; }

    void updateMesh(const std::vector<BasicPoint> & vertices , const std::vector<Tetrahedron> & tetrahedra, float xshearstep );
    void updateVertices( std::vector<BasicPoint> & points );
    void computePseudoRasterDeformation(VoxelGrid & HRGrid, VoxelGrid & HRShearedGrid, QString & foldername, VoxelGrid & HRresult, const QString & outputFoldername, int threshold);
    void computePseudoRasterDeformation(VoxelGrid & HRGrid, VoxelGrid & HRShearedGrid, std::ifstream & inputImaFile, VoxelGrid & HRresult, const QString & outputFoldername, int threshold, int interpolation, QStatusBar *statusbar );
    void computeBarycentricCoordinates(int tet_id , const BasicPoint & point, float & ld0, float & ld1, float & ld2, float & ld3 );
    void computePseudoRasterDeformation(VoxelGrid & HRGrid, VoxelGrid & HRShearedGrid,std::ifstream & inputImaFile, VoxelGrid & HRresult, std::ofstream & outputImaFile, int threshold );
    int getValueFromFile( int i, int j, int k, VoxelGrid & HRGrid, VoxelGrid & HRShearedGrid, std::ifstream & inputImaFile, float maxshear, float z_shear);
protected :
    void computeBB();

    std::vector<BasicPoint> _vertices;
    std::vector<Tetrahedron> _tetrahedra;

    BasicPoint BBMin;
    BasicPoint BBMax;
    bool dataLoaded;
    std::vector<BasicPoint> initial_vertices_positions;
    VoxelGrid * voxel_grid;

    float _xshearstep;

};

#endif // TETMESHCREATOR_H
