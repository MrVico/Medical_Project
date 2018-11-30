#ifndef VOXELGRID_H
#define VOXELGRID_H

#include <map>
#include <vector>
#include <iostream>
#include "BasicPoint.h"
#include "Voxel.h"
#include <QColor>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define DEFAULT_VOFFSET 0

typedef int Subdomain_index ;
typedef char GRID_TYPE;

// -------------------------------------------------
// Intermediate Voxel structure for hashed adjacency
// -------------------------------------------------

class VoxelGrid
{
public:

    VoxelGrid() : _offset(BasicPoint(0.,0.,0.))
    {
        init(0, 0, 0, 1.,1.,1.);
    }

    // ~VoxelGrid(){clear(); }
    //TODO destructor!!!

    VoxelGrid(const std::vector<GRID_TYPE> & data, unsigned int x,unsigned int y, unsigned int z, float _dx, float _dy, float _dz,
              BasicPoint offset = BasicPoint(0.,0.,0.), int _VOffset = DEFAULT_VOFFSET);

    VoxelGrid(const std::vector<GRID_TYPE> & data, unsigned int x,unsigned int y, unsigned int z, float _dx, float _dy, float _dz,bool sort,
              BasicPoint offset = BasicPoint(0.,0.,0.), int _VOffset = DEFAULT_VOFFSET);


    VoxelGrid(unsigned int nx,unsigned int ny, unsigned int nz, float dx, float dy, float dz, bool allocate_memory, BasicPoint offset = BasicPoint(0.,0.,0.), int _VOffset = DEFAULT_VOFFSET):
        _offset(offset),
        VOffset(_VOffset){
        init(nx, ny, nz, dx, dy, dz, allocate_memory);
    }

    VoxelGrid(const VoxelGrid & voxelGrid){
        init(voxelGrid._dim[0], voxelGrid._dim[1], voxelGrid._dim[2], voxelGrid._d[0], voxelGrid._d[1], voxelGrid._d[2]);

        voxels = voxelGrid.voxels;
        positions = std::vector<BasicPoint>(voxelGrid.positions);
        sortedVoxels = voxelGrid.sortedVoxels;
        subdomain_indices = voxelGrid.subdomain_indices;
        colors = voxelGrid.colors;
        BBMin = voxelGrid.BBMin;
        BBMax = voxelGrid.BBMax;
        VOffset = voxelGrid.VOffset;
        clippingNormal = voxelGrid.clippingNormal;
        pointOnClipping = voxelGrid.pointOnClipping;

        _max = voxelGrid._max;
        _size = voxelGrid._size;

        grid_data = voxelGrid.grid_data;

        vnumbers = voxelGrid.vnumbers;
    }

    inline float square( float x ){ return x*x ;}

    inline float wendland_kernel( float x , float h )
    {
        if( x > h )
            return 0.f;

        return square( square( 1.f - x/h))*(4*x/h+1);
    }


    void computeEnvelopVoxels();
    int getIndice(int i, int j , int k);

    unsigned int getSubdomainNumber(){ return subdomain_indices.size(); }
    void getSubdomainIndices(std::vector<Subdomain_index> & _subdomain_indices);
    GRID_TYPE getGridTypeValue(Subdomain_index si);

    Subdomain_index getSubdomainIndex(GRID_TYPE g_value);

    std::vector< Voxel > & getVoxels(){ return voxels ; }
    std::vector< Voxel > const & getVoxels() const { return voxels ; }

    std::vector< BasicPoint > & getPositions(){ return positions ; }
    std::vector< BasicPoint > const & getPositions() const { return positions ; }

    const BasicPoint & getOffset() const { return _offset; }

    bool is_in_domains(unsigned int i, unsigned int j, unsigned int k){ GRID_TYPE v = value(i,j,k); if( v <= BACKGROUND_GRID_VALUE ) return false; return true; }
    bool is_in_offseted_domains(unsigned int i, unsigned int j, unsigned int k){ GRID_TYPE v = value(i,j,k); if( v < BACKGROUND_GRID_VALUE ) return false; return true; }
    bool is_in_offseted_domains(unsigned int i){ GRID_TYPE v = value(i); if( v < BACKGROUND_GRID_VALUE ) return false; return true; }
    unsigned int xdim(){ return _dim[0]; }
    unsigned int ydim(){ return _dim[1]; }
    unsigned int zdim(){ return _dim[2]; }

    unsigned int size(){ return _size; }

    unsigned int dim(int i){ assert(i >= 0 && i < 3); return _dim[i];}

    double dx(){ return _d[0]; }
    double dy(){ return _d[1]; }
    double dz(){ return _d[2]; }

    const BasicPoint & offset(){ return _offset ; }

    std::vector<GRID_TYPE> & data(){ return grid_data; }
    const std::vector<GRID_TYPE> & data() const { return grid_data; }

    double d(int i){ assert(i >= 0 && i < 3); return _d[i];}

    void getDistanceFiels( std::vector<float> &d_field );
    float getRadius(){return (BBMax-BBMin).norm()/2.;}
    BasicPoint getCenter(){return (BBMax-BBMin)/2.;}

    const BasicPoint & getBBMax(){return BBMax;}
    const BasicPoint & getBBMin(){return BBMin;}
    const BasicPoint & getOffset(){return _offset;}
    int getVOffset(){return VOffset; }
    const Voxel & getVoxelOffset(){return voffset; }
    void computeBoundingBox();
    void setClipeEquation(const BasicPoint & clipN, const BasicPoint & pointN){ clippingNormal = clipN; clippingNormal.normalize(); pointOnClipping = pointN; }
    void setOrthongonalCut(const BasicPoint & _cut, const BasicPoint & _cutDirection){ cut = _cut; cutDirection = _cutDirection; }
    void rasterize(VoxelGrid & deformed_grid, double resolution);
    void rasterizeClosest(VoxelGrid & deformed_grid, unsigned int neighborhoodSize, double resolution);
    void rasterizeUsingVectorFields(VoxelGrid & deformed_grid, unsigned int neighborhoodSize, double resolution);

    void dilatation( int voffsetNb , unsigned int dx, unsigned int dy, unsigned int dz, const Voxel & voffset, bool s_changed);

    void downSample();
    GRID_TYPE value( unsigned int i ) const ;
    GRID_TYPE originalValue( unsigned int i ) { return value(i) - BACKGROUND_GRID_VALUE; }

    void setSegmentation(const std::vector<GRID_TYPE> &segmentation);
    void saveInfoToFile(const std::string & filename);

    int getVoxelNumber(const Subdomain_index si );

    void compareVolumeLossWithGrid(  VoxelGrid & grid , std::vector< float > & weight_diff);
    void computeElongation(const std::vector<BasicPoint> & defPositions, std::vector<int> & values, float & min_value, float & max_value);

    void clear();
    void clearAll(){
        grid_data.clear();
        clear();
    }

    void setColors( const std::vector<QColor> & _colors );
    void setWhite();

    void computeDefaultColors();

    void checkForOutliers();
    bool ignoreOutliers();
    void drawOutliers();
    void computeDistortion(const std::vector<BasicPoint> & defPositions, std::vector<int> & values, float & min_value, float & max_value);
    void computeDistortion(const std::vector<BasicPoint> &defPositions, std::vector<float> &distortions);
    void getDeformedVoxelGrid(VoxelGrid & voxelGrid, const BasicPoint & bbmin, const BasicPoint & bbmax);

    void projectInRegularGrid(VoxelGrid & deformed_grid, std::map<unsigned int, unsigned int> & VMap, double resolution=1., bool sort = false);

    float computeMinPointDist();

    std::vector <QColor> & getColors(){return colors; }
    BasicPoint getWorldCoordinate( int i, int j, int k );
    BasicPoint getWorldCoordinateVBottom( unsigned int i, unsigned int j, unsigned int k );
    void setValue( unsigned int i, unsigned int j, unsigned int k, const GRID_TYPE & v );
    void sort();
    void sortSeg();

    GRID_TYPE value(const BasicPoint & point )const { Voxel v = getGridCoordinate(point); if(isInGrid(v))return value(v); return DEFAULT_GRID_VALUE;}
    GRID_TYPE value( unsigned int i, unsigned int j, unsigned int k ) const;

    void computePositions();
    unsigned int getGridIndex(unsigned int i, unsigned int j, unsigned int k ) const;
    unsigned int getZFirstGridIndex(unsigned int i, unsigned j, unsigned int k ) const;
    void clearToDefault();
    Voxel getGridCoordinate(const BasicPoint & point) const ;
    Voxel getGridCoordinateNN(const BasicPoint & point) const;
    void getNeighborhoodBorns(const Voxel & voxel, Voxel & v_born_min, Voxel & v_born_max, int kernel_size = 1 );
    unsigned int getGridIndex(const Voxel & voxel ) const ;
    unsigned int getZFirstGridIndex(const Voxel & voxel ) const ;
    void addAndUpdateVoxel(Voxel & v , GRID_TYPE i){sortedVoxels[i].push_back(voxels.size()); voxels.push_back(v); positions.push_back(getWorldCoordinate(v)); }
    void setValue( unsigned int i, const GRID_TYPE & v );
    bool isInGrid(const Voxel & voxel)const;

    Subdomain_index getInitialSubdomainIndex(unsigned int i){
        return value(i) - BACKGROUND_GRID_VALUE;
    }

    bool updateOffset( int _Voffset, VoxelGrid & segmentation );
    static GRID_TYPE DEFAULT_GRID_VALUE;
    static GRID_TYPE BACKGROUND_GRID_VALUE;
    void sortAndInit(const std::vector<GRID_TYPE> & data);

    void removeUnconnectedVoxels();
    void keepLabelsValues(  std::map<Subdomain_index, bool> & label_to_keep, VoxelGrid & grid );
    void separateInConnectedComponents( std::map<Subdomain_index, bool> & label_to_keep, VoxelGrid & grid );
    void dilateToRemoveVisible( std::map<Subdomain_index, bool> & label_to_remove, VoxelGrid & grid );
    void keepVisiblePart( VoxelGrid & grid);

    bool changeLabelValue( Subdomain_index from, Subdomain_index to );
    void computeVNumbers( std::map<Subdomain_index, int> & vNbs);
    int counter;

    std::vector<Voxel> & getEnvelopVoxels(){return envelopVoxels; }
    void getSegmentation(VoxelGrid & segmentation);
    void automaticSeparation();
    bool isPairToSeparate( Subdomain_index _l1, Subdomain_index _l2 );
    void replaceVisibleLabelsValue( std::map<Subdomain_index, bool> & labels_to_decompose, Subdomain_index to, VoxelGrid & grid );
    void fitToData( int offset , VoxelGrid & grid );
    void addLayer( int id , VoxelGrid & grid );

    void setColors(    std::map<Subdomain_index, QColor> _fileColors) {fileColors =_fileColors;}

    BasicPoint getWorldCoordinate(const Voxel & voxel);
    BasicPoint getWorldCoordinateVBottom(const Voxel & voxel);
protected:

    std::vector<unsigned int> offsets;
    std::vector<Voxel> envelopVoxels;
    std::vector<std::pair<int, int> > VBOInfos;
    std::vector<int> userOffsets;

    std::map<Subdomain_index, QColor> fileColors;

    void init( unsigned int nx, unsigned int ny, unsigned int nz, float dx, float dy, float dz, bool allocate_memory = true );
    void clearToZero();



    bool ajustSizeToOffset( );

    void getDeformedVoxelGrid(VoxelGrid & voxelGrid, double resolution);


    bool dilatate();
    bool dilatate(const Voxel & voxel);
    bool dilatate(unsigned int i , unsigned int j, unsigned int k);
    bool dilatation( const std::vector<int> & envelopVoxels );

    GRID_TYPE value( const Voxel & voxel ) const ;

    void setValue( const Voxel & voxel, const GRID_TYPE & v );



    void getNeighborhoodBorns(int i, int j, int k, Voxel & v_born_min, Voxel & v_born_max, int kernel_size = 1 );

    void collectSixNeighborhoodVoxels(const Voxel & voxel, std::vector<Voxel> & neighbors);
    void collectSixNeighborhoodVoxels(unsigned int i, unsigned int j, unsigned int k, std::vector<Voxel> & neighbors);


    bool isInGrid(int i, int j, int k)const;
    bool isInGrid(int value)const;



    BasicPoint cut;
    BasicPoint cutDirection;

    void drawCube(unsigned int vId);
    void drawCube(const BasicPoint & center);
    bool isVisiblePoint( const BasicPoint & point );

    void findGoodTexSize( unsigned int uiNumValuesNeeded, int& iTexWidth, int& iTexHeight, bool bUseAlphaChannel );

    int getNumberOfBoundaryVoxels();

    int round(float r);

    std::vector<GRID_TYPE> grid_data;
    std::vector<Voxel> voxels;
    std::vector<BasicPoint> positions;

    std::vector<std::vector<int> > sortedVoxels;

    std::vector<GRID_TYPE> subdomain_indices;
    std::vector<Subdomain_index> initial_subdomain_indices;
    std::vector<int> vnumbers;
    std::vector<QColor> colors;

    BasicPoint BBMin;
    BasicPoint BBMax;

    Voxel Vmin;
    Voxel Vmax;

    BasicPoint clippingNormal;
    BasicPoint pointOnClipping;

    BasicPoint _offset;

    unsigned int _dim[3];
    float _d[3];
    unsigned int _size;
    float _max;

    int VOffset;

    Voxel voffset;

    std::vector<int> outliers;

};

#endif // VOXELGRID_H
