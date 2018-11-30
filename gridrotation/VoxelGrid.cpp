#include "VoxelGrid.h"
//#include "GLUtilityMethods.h"
#include "float.h"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <queue>
#include <QTime>
#include <cmath>
#include <climits>

GRID_TYPE VoxelGrid::DEFAULT_GRID_VALUE = -125;
GRID_TYPE VoxelGrid::BACKGROUND_GRID_VALUE = -124;

VoxelGrid::VoxelGrid(const std::vector<GRID_TYPE> & data, unsigned int nx, unsigned int ny, unsigned int nz, float dx, float dy, float dz, BasicPoint offset, int _VOffset):
    _offset(offset),
    VOffset(_VOffset)

{
    init(nx, ny, nz, dx, dy, dz);
    
    sortAndInit(data);
    
}

VoxelGrid::VoxelGrid(const std::vector<GRID_TYPE> & data, unsigned int nx, unsigned int ny, unsigned int nz, float dx, float dy, float dz, bool sorting, BasicPoint offset, int _VOffset):
    _offset(offset),
    VOffset(_VOffset)
{
    init(nx, ny, nz, dx, dy, dz);
    
    grid_data = data;
    if(sorting)
        this->sort();
}

void VoxelGrid::init(unsigned int nx,unsigned int ny, unsigned int nz, float dx, float dy, float dz, bool allocate_memory){
    
    cut = BasicPoint(0.,0.,0.);
    cutDirection = BasicPoint(1.,1.,1.);

    // counter = 50;
    _dim[0] = nx;
    _dim[1] = ny;
    _dim[2] = nz;
    
    _size = nx*ny*nz;
    
    _d[0] = dx;
    _d[1] = dy;
    _d[2] = dz;
    
    _max = std::max(nx*dx, ny*dy);
    _max = std::max(nz*dz, _max);
    
    if(allocate_memory){
        grid_data.clear();
        grid_data.resize(_size,DEFAULT_GRID_VALUE);
    }

}

void VoxelGrid::clearToDefault(){
    grid_data.clear();
    grid_data.resize(_size, DEFAULT_GRID_VALUE);
}

void VoxelGrid::clearToZero(){
    grid_data.clear();
    grid_data.resize(_size, BACKGROUND_GRID_VALUE);
}

void VoxelGrid::downSample( ){
    unsigned int n[3];
    unsigned int p[3];
    for(int i = 0 ; i < 3 ; i ++){
        n[i] = int(_dim[i]/2.);
        p[i]=0;
        if( n[i]%2 > 0) p[i]++;
    }
    
    
    VoxelGrid grid = VoxelGrid(n[0]+p[0], n[1]+p[1], n[2]+p[2], _d[0]*2, _d[1]*2, _d[2]*2, true, _offset );
    grid.clearToZero();
    
    for( unsigned int i = 0, oi = 0 ; i < n[0] ; i++, oi+=2 ){
        for( unsigned int j = 0, oj = 0 ; j < n[1] ; j++, oj+=2){
            for( unsigned int k = 0, ok = 0 ; k < n[2] ; k++, ok += 2 ){
                
                std::map<int, int> ids;
                
                for( unsigned int x = 0; x < 2 ; x++ ){
                    for( unsigned int y = 0; y < 2 ; y++ ){
                        for( unsigned int z = 0; z < 2 ; z++ ){
                            ids[std::max(grid_data[getGridIndex(oi+x,oj+y,ok+z)], BACKGROUND_GRID_VALUE)]++;
                        }
                    }
                }
                
                int maxNb = 0 ;
                int id;
                
                for(std::map<int, int>::iterator it = ids.begin() ; it != ids.end() ; it++){
                    if( it->second > maxNb ){
                        id = it->first;
                        maxNb = it->second;
                    }
                }
                
                grid.setValue( i, j, k, id );
            }
        }
    }
    
    clearAll();
    
    _offset = grid.getOffset();
    init(grid.dim(0), grid.dim(1), grid.dim(2), grid.d(0), grid.d(1), grid.d(2));
    
    sortAndInit(grid.data());
    
}


void VoxelGrid::computeVNumbers( std::map<Subdomain_index, int> & vNbs){
    
    vNbs.clear();
    for(unsigned int i = 0 ; i < _size ; i++){
        
        GRID_TYPE v_value = value(i);
        if( v_value == DEFAULT_GRID_VALUE ){
            v_value = BACKGROUND_GRID_VALUE;
        }
        
        Subdomain_index si = getSubdomainIndex(v_value);
        if(si > 0)
            vNbs[si] ++;
    }
}


void VoxelGrid::sortSeg(){

    voxels.clear();
    positions.clear();
    sortedVoxels.clear();
    subdomain_indices.clear();
    colors.clear();
    vnumbers.clear();

    for(unsigned int i = 0 ; i < _dim[0] ; i++){
        for(unsigned int j = 0 ; j < _dim[1] ; j++){
            for(unsigned int k = 0 ; k < _dim[2] ; k++){

                int index = getGridIndex(i,j,k);

                GRID_TYPE v_value = value(index);
                if( v_value == DEFAULT_GRID_VALUE ){
                    setValue( index, BACKGROUND_GRID_VALUE );
                    v_value = BACKGROUND_GRID_VALUE;
                }

                if(v_value> BACKGROUND_GRID_VALUE){
                    std::vector<GRID_TYPE>::iterator it = std::find(subdomain_indices.begin(), subdomain_indices.end(), v_value);

                    if( it == subdomain_indices.end() ){
                        subdomain_indices.push_back(v_value);
                        vnumbers.push_back(1);
                        sortedVoxels.push_back(std::vector<int>(1, voxels.size()));

                    } else {
                        int id = it - subdomain_indices.begin();
                        sortedVoxels[id].push_back(voxels.size());
                        vnumbers[id] +=1;
                    }
                    Voxel voxel = Voxel(i,j,k);
                    voxels.push_back(voxel);

                    Vmin = Voxel::min(Vmin, voxel);
                    Vmax = Voxel::max(Vmax, voxel);

                    std::vector<Voxel> neighbors;
                    collectSixNeighborhoodVoxels(i, j, k, neighbors);

                }
            }
        }
    }

    //   computeDefaultColors();

    //computePositions();

    //createDisplayLists();

    BasicPoint max (dx()/2., dy()/2., dz()/2.);
    BasicPoint min (-dx()/2., -dy()/2., -dz()/2.);

    //    initTextures( min, max );
    //    computeTranslationData();
    std::cout << "Sorted " << std::endl;
}

void VoxelGrid::sort(){
    
    voxels.clear();
    positions.clear();
    sortedVoxels.clear();
    subdomain_indices.clear();
    colors.clear();
    vnumbers.clear();
    
    for(unsigned int i = 0 ; i < _dim[0] ; i++){
        for(unsigned int j = 0 ; j < _dim[1] ; j++){
            for(unsigned int k = 0 ; k < _dim[2] ; k++){
                
                int index = getGridIndex(i,j,k);
                
                GRID_TYPE v_value = value(index);
                if( v_value == DEFAULT_GRID_VALUE ){
                    setValue( index, BACKGROUND_GRID_VALUE );
                    v_value = BACKGROUND_GRID_VALUE;
                }
                
                if(v_value> BACKGROUND_GRID_VALUE){
                    std::vector<GRID_TYPE>::iterator it = std::find(subdomain_indices.begin(), subdomain_indices.end(), v_value);

                    if( it == subdomain_indices.end() ){
                        subdomain_indices.push_back(v_value);
                        vnumbers.push_back(1);
                        sortedVoxels.push_back(std::vector<int>(1, voxels.size()));

                    } else {
                        int id = it - subdomain_indices.begin();
                        sortedVoxels[id].push_back(voxels.size());
                        vnumbers[id] +=1;
                    }

                    voxels.push_back(Voxel(i,j,k));
                }
            }
        }
    }

    computeDefaultColors();

    computePositions();

    std::cout << "Sorted " << std::endl;
}

bool VoxelGrid::updateOffset( int _Voffset, VoxelGrid & segmentation ){
    
    if (_Voffset == VOffset) return false;
    
    std::cout << "Updating voxel offset " << std::endl;
    bool s_change = false;
    if( _Voffset > VOffset ){
        VOffset = _Voffset;
        if( ajustSizeToOffset() ){
            //std::cout << _Voffset - VOffset << std::endl;
            computePositions();
            computeBoundingBox();
            s_change = true;
        }
    }
    
    segmentation.dilatation( _Voffset, _dim[0], _dim[1], _dim[2],voffset, false );
    
    return s_change;
#if 0
    
    if (_Voffset == VOffset) return;
    
    std::cout << "Updating voxel offset " << std::endl;
    
    if( _Voffset < 0 ){
        std::cout << _Voffset - VOffset << std::endl;
        clear();
        dilatation(_Voffset);
        VOffset = 1;
        sortAndInit(grid_data);
        VOffset = _Voffset;
        
    }
    
    else {
        
        std::vector<GRID_TYPE>::iterator it = std::find(subdomain_indices.begin(), subdomain_indices.end(), BACKGROUND_GRID_VALUE);
        
        int id = 0;
        if( it == subdomain_indices.end() ){
            
            subdomain_indices.push_back(BACKGROUND_GRID_VALUE);
            id = sortedVoxels.size();
            sortedVoxels.push_back(std::vector<int> ());
            
        } else {
            
            id = it - subdomain_indices.begin();
            
            // std::cout <<"Supposed to be zero " << (int)value(voxels[sortedVoxels[id].front()]) << std::endl;
            for(unsigned int i = 0 ; i < sortedVoxels[id].size() ; i++){
                setValue(voxels[sortedVoxels[id][i]], DEFAULT_GRID_VALUE);
            }
            
            
            voxels.erase( voxels.begin()+sortedVoxels[id].front(), voxels.begin()+sortedVoxels[id].front() + sortedVoxels[id].size() );
            positions.erase( positions.begin()+sortedVoxels[id].front(), positions.begin()+sortedVoxels[id].front() + sortedVoxels[id].size() );
            sortedVoxels[id].clear();
            
        }
        
        VOffset = _Voffset;
        std::cout << "Performing dilatation " << std::endl;
        if( dilatation(envelopVoxels) ){
            
            computePositions();
            clearDisplayLists();
            createDisplayLists();
            
        } else {
            
            for(unsigned int i = 0 ; i < sortedVoxels[id].size() ; i++){
                const Voxel & voxel = voxels[sortedVoxels[id][i]];
                positions.push_back(getWorldCoordinate(voxel));
            }
            glDeleteLists(displayListIds[id], 1);
            createDisplayList(id);
        }
#endif
        
    }
    
    
    void VoxelGrid::sortAndInit(const std::vector<GRID_TYPE> & data){
        
        clear();
        envelopVoxels.clear();
        
        unsigned int nx = xdim();
        unsigned int ny = ydim();
        unsigned int nz = zdim();
        
        Vmin = Voxel(nx, ny, nz);
        Vmax = Voxel(0, 0, 0);
        /*
        std::vector<GRID_TYPE> data ( tmp.size() );
        for(unsigned int i = 0, inv_i = nx-1 ; i < nx ; i++, inv_i--){
            for(unsigned int j = 0 ; j < ny ; j++){
                for(unsigned int k = 0 ; k < nz ; k++){
                    data[ getGridIndex(inv_i,j,k) ] = tmp[ getGridIndex(i,j,k) ];
                }
            }
        }
*/
        //  subdomain_indices.push_back(BACKGROUND_GRID_VALUE);
        //  sortedVoxels.push_back(std::vector<int>());
        for(unsigned int i = 0 ; i < nx ; i++){
            for(unsigned int j = 0 ; j < ny ; j++){
                for(unsigned int k = 0 ; k < nz ; k++){
                    
                    int indice = getGridIndex(i, j, k);
                    GRID_TYPE si = data[indice];
                    
                    int current_i = voxels.size();
                    Voxel voxel(i, j, k);
                    
                    if( si > BACKGROUND_GRID_VALUE ){
                        
                        std::vector<GRID_TYPE>::iterator it = std::find(subdomain_indices.begin(), subdomain_indices.end(), si);
                        int id = it - subdomain_indices.begin();

                        if( it == subdomain_indices.end() ){
                            subdomain_indices.push_back(si);
                            sortedVoxels.push_back(std::vector<int> (1,voxels.size()));
                        } else {
                            sortedVoxels[id].push_back(voxels.size());
                        }
                        
                        voxels.push_back(voxel);
                        
                        Vmin = Voxel::min(Vmin, voxel);
                        Vmax = Voxel::max(Vmax, voxel);
                        
                        std::vector<Voxel> neighbors;
                        collectSixNeighborhoodVoxels(i, j, k, neighbors);
                        
                        for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                            const Voxel & n_voxel = neighbors[n];
                            
                            if( data[getGridIndex( n_voxel )] <= BACKGROUND_GRID_VALUE || i ==0 || j ==0 || k ==0 || i == nx-1 || j == ny-1 || k == nz-1   ){
                                envelopVoxels.push_back(Voxel(i,j,k));
                                break;
                            }
                        }
                        
                    }
                    
                    if(si == BACKGROUND_GRID_VALUE)
                        si = DEFAULT_GRID_VALUE;
                    
                    grid_data[indice] = si;
                    
                }
            }
        }
        std::cout << "sorted " << std::endl;
        //    dilatation(envelopVoxels);
        
        
        for(unsigned int i = 0 ; i < sortedVoxels.size() ; i++){
            vnumbers.push_back(sortedVoxels.size());
            //        std::cout << "si : " << (int)subdomain_indices[i] << std::endl;
        }
        
        
        
        std::cout << "Dilatation done, updating positions " << std::endl;
        computeDefaultColors();
        
        //_offset -= getWorldCoordinate(Vmin);
        computePositions();
        std::cout << "Grid built " << positions.size() << std::endl;
        
        //    initTextures( min, max );
        //    computeTranslationData();
        
    }
    
    int VoxelGrid::getNumberOfBoundaryVoxels(){
        
        int count = 0 ;
        for( unsigned int i = 0 ; i < voxels.size() ; i++ ){
            const Voxel & voxel = voxels[i];
            std::vector<Voxel> neighbors;
            GRID_TYPE gd = value(voxel);
            collectSixNeighborhoodVoxels(voxel, neighbors);
            
            for( unsigned int j = 0 ; j < neighbors.size() ; j++ ){
                if( gd != value(neighbors[j]) ){
                    count ++;
                    break;
                }
            }
            
        }
        
        return count;
    }
    
    void VoxelGrid::saveInfoToFile(const std::string &filename){
        std::ofstream myfile;
        myfile.open(filename.c_str());
        if (!myfile.is_open())
        {
            std::cout << filename << " cannot be opened" << std::endl;
            return;
        }
        
        myfile << "Dimensions" <<std::endl;
        
        myfile << _dim[0] << " " << _dim[1] << " " << _dim[2] << std::endl;
        myfile << _d[0] << " " << _d[1] << " " << _d[2] << std::endl;
        
        
        myfile << "Size on memory of " << voxels.size() << " voxels " << std::endl;
        int globalDataSize = 0;
        int tmp = sizeof(std::vector<GRID_TYPE>) + grid_data.capacity() * sizeof(GRID_TYPE);
        globalDataSize += tmp ;
        myfile << "Grid data " << tmp << std::endl;
        tmp = sizeof(std::vector<Voxel>) + voxels.capacity() * sizeof(Voxel);
        globalDataSize += tmp ;
        myfile << "Voxels " << tmp << std::endl;
        tmp = sizeof(std::vector<BasicPoint>) + positions.capacity() * sizeof(BasicPoint);
        globalDataSize += tmp ;
        myfile << "Positions " << tmp << std::endl;
        tmp = sizeof(std::vector<int>) + voxels.capacity() * sizeof(int);
        globalDataSize += tmp ;
        myfile << "Sorted voxel " << tmp << std::endl;
        myfile << "Total : " << globalDataSize << std::endl;
        
        
        myfile << std::endl;
        int nbBoundaryVoxels = getNumberOfBoundaryVoxels();
        myfile << nbBoundaryVoxels << " voxels on boundaries"<< std::endl ;
        myfile << float(nbBoundaryVoxels)*100./voxels.size() << " % "<< std::endl ;
        myfile << std::endl;
        myfile.close ();
        
    }
    
    void VoxelGrid::findGoodTexSize( unsigned int uiNumValuesNeeded, int& iTexWidth, int& iTexHeight, bool bUseAlphaChannel )
    {
        unsigned int uiFactor = ( bUseAlphaChannel ? 4 : 3 );
        int iRoot = ( int )sqrtf( ( float ) ( uiNumValuesNeeded / uiFactor ) );
        
        iTexWidth = iRoot + 1;
        iTexHeight = iRoot;
        if ( ( unsigned int ) ( iTexWidth * iTexHeight ) * uiFactor < uiNumValuesNeeded ) {
            ++iTexHeight;
        }
        if ( ( unsigned int ) ( iTexWidth * iTexHeight ) * uiFactor < uiNumValuesNeeded ) {
            ++iTexWidth;
        }
    }
    


    void VoxelGrid::clear(){
        
        voxels.clear();
        positions.clear();
        sortedVoxels.clear();
        subdomain_indices.clear();
        colors.clear();
        vnumbers.clear();
        envelopVoxels.clear();
        
        
    }

    
    bool VoxelGrid::dilatate(){
        bool result = false;
        const std::vector<GRID_TYPE> currentState (grid_data);
        for(unsigned int i = 0 ; i < _dim[0] ; i++){
            for(unsigned int j = 0 ; j < _dim[1] ; j++){
                for(unsigned int k = 0 ; k < _dim[2] ; k++){
                    int label = currentState[getGridIndex(i,j,k)];
                    if( label >= BACKGROUND_GRID_VALUE ){
                        bool stepResult = dilatate(i, j, k);
                        if( stepResult && label > BACKGROUND_GRID_VALUE )
                            result = true;
                    }
                }
            }
        }
        return result;
    }
    
    
    bool VoxelGrid::dilatate(const Voxel & voxel){
        return dilatate(voxel.i(), voxel.j(), voxel.k());
    }
    
    bool VoxelGrid::dilatate(unsigned int i , unsigned int j, unsigned int k){
        
        GRID_TYPE label = value(i, j, k);
        if( label <= BACKGROUND_GRID_VALUE ) return false;
        
        bool result = false;
        std::vector<Voxel> neighbors;
        collectSixNeighborhoodVoxels(i, j, k, neighbors);
        
        for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
            const Voxel & voxel = neighbors[n];
            
            if( value( voxel ) <= BACKGROUND_GRID_VALUE ){
                setValue(voxel, label);
                result = true;
            }
        }
        
        return result;
    }
    
    unsigned int VoxelGrid::getGridIndex(const Voxel & voxel ) const {
        return getGridIndex(voxel.i(), voxel.j(), voxel.k());
    }
    
    unsigned int VoxelGrid::getGridIndex(unsigned int i, unsigned j, unsigned int k ) const {
        return i + j*_dim[0] + k*_dim[0]*_dim[1];
    }
    
    unsigned int VoxelGrid::getZFirstGridIndex(const Voxel & voxel ) const {
        return getZFirstGridIndex(voxel.i(), voxel.j(), voxel.k());
    }

    unsigned int VoxelGrid::getZFirstGridIndex(unsigned int i, unsigned j, unsigned int k ) const {
        return k + i*_dim[2] + j*_dim[2]*_dim[0];
    }

    GRID_TYPE VoxelGrid::value( const Voxel & voxel )const {
        return value(voxel.i(), voxel.j(), voxel.k());
    }
    
    GRID_TYPE VoxelGrid::value( unsigned int i, unsigned int j, unsigned int k ) const {
        return value(getGridIndex(i, j, k));
    }
    
    GRID_TYPE VoxelGrid::value( unsigned int i ) const {
        return grid_data[i];
    }
    
    void VoxelGrid::setValue( const Voxel & voxel, const GRID_TYPE & v ){
        setValue( voxel.i(), voxel.j(), voxel.k(), v );
    }
    
    void VoxelGrid::setValue( unsigned int i, unsigned int j, unsigned int k, const GRID_TYPE & v ){
        setValue(getGridIndex(i, j, k), v);
    }
    
    void VoxelGrid::setValue( unsigned int i, const GRID_TYPE & v ){
        grid_data[i] = v;
    }
    
    void VoxelGrid::getNeighborhoodBorns(const Voxel & voxel, Voxel & v_born_min, Voxel & v_born_max, int kernel_size ){
        getNeighborhoodBorns(voxel.i(), voxel.j(), voxel.k(), v_born_min, v_born_max, kernel_size);
    }
    
    void VoxelGrid::getNeighborhoodBorns(int i, int j, int k, Voxel & v_born_min, Voxel & v_born_max, int kernel_size ){
        v_born_min = Voxel::max( Voxel( i-kernel_size, j-kernel_size, k-kernel_size ), Voxel( 0,0,0 ) );
        v_born_max = Voxel::min( Voxel( i+kernel_size, j+kernel_size, k+kernel_size ), Voxel( _dim[0]-1, _dim[1]-1, _dim[2]-1 ) );
    }
    
    void VoxelGrid::collectSixNeighborhoodVoxels(const Voxel & voxel, std::vector<Voxel> & neighbors){
        collectSixNeighborhoodVoxels(voxel.i(), voxel.j(), voxel.k(), neighbors);
    }
    
    void VoxelGrid::collectSixNeighborhoodVoxels(unsigned int i, unsigned int j, unsigned int k, std::vector<Voxel> & neighbors){
        
        if((int)i-1 >= 0)
            neighbors.push_back( Voxel((int)i-1, j, k) );
        if(i+1 < _dim[0])
            neighbors.push_back( Voxel(i+1, j, k) );
        
        if((int)j-1 >= 0)
            neighbors.push_back( Voxel(i, ((int)j-1), k) );
        if(j+1 < _dim[1])
            neighbors.push_back( Voxel(i, (j+1), k) );
        
        if((int)k-1 >= 0)
            neighbors.push_back( Voxel(i, j, ((int)k-1)) );
        if(k+1 < _dim[2])
            neighbors.push_back( Voxel(i, j, (k+1)) );
    }
    
    bool VoxelGrid::isInGrid(const Voxel & voxel)const {
        return isInGrid(voxel.i(), voxel.j(), voxel.k());
    }
    
    bool VoxelGrid::isInGrid(int i, int j, int k)const{
        if( i < 0 || i >= _dim[0] || j < 0 || j >= _dim[1] ||k >= _dim[2] ||k < 0 )
            return false;
        return isInGrid(getGridIndex(i,j,k));
    }
    
    bool VoxelGrid::isInGrid(int value)const {
        return (value >= 0) && (value < grid_data.size());
    }
    
    BasicPoint VoxelGrid::getWorldCoordinateVBottom(const Voxel & voxel){
        return getWorldCoordinateVBottom( voxel.i(), voxel.j(), voxel.k() ) ;
    }

    BasicPoint VoxelGrid::getWorldCoordinateVBottom( unsigned int i, unsigned int j, unsigned int k ){
        // return BasicPoint((i+0.5 -voffset.i())*_d[0], (j+0.5-voffset.j())*_d[1], (k+0.5-voffset.k())*_d[2])+_offset ;
        return BasicPoint((i+0.5)*_d[0], (j+0.5)*_d[1], k*_d[2])+_offset ;
    }

    BasicPoint VoxelGrid::getWorldCoordinate(const Voxel & voxel){
        return getWorldCoordinate( voxel.i(), voxel.j(), voxel.k() ) ;
    }
    
    BasicPoint VoxelGrid::getWorldCoordinate( int i, int j, int k ){
        // return BasicPoint((i+0.5 -voffset.i())*_d[0], (j+0.5-voffset.j())*_d[1], (k+0.5-voffset.k())*_d[2])+_offset ;
        return BasicPoint((i+0.5)*_d[0], (j+0.5)*_d[1], (k+0.5)*_d[2])+_offset ;
    }

    Voxel VoxelGrid::VoxelGrid::getGridCoordinate(const BasicPoint & point) const {
        BasicPoint pos = point - _offset;
        return Voxel(int(pos[0]/_d[0]) , int(pos[1]/_d[1]), int(pos[2]/_d[2]));
    }
    
    Voxel VoxelGrid::VoxelGrid::getGridCoordinateNN(const BasicPoint & point) const {
        BasicPoint pos = point - _offset;
        pos = BasicPoint(pos[0]/_d[0] , pos[1]/_d[1], pos[2]/_d[2]);
        return Voxel(int(pos[0]) , int(pos[1]), int(pos[2]+0.5*_d[2]));
    }

    int VoxelGrid::round(float r){
        return (r - std::floor(r) > 0.5) ? std::ceil(r) : std::floor(r);
    }
    
    void VoxelGrid::computeDefaultColors(){

        colors.clear();
        colors.resize(subdomain_indices.size(), QColor(0,0,0) );
        int off = 0;
        int max_si = 0;
        for( unsigned int i = 0 ; i < subdomain_indices.size() ; i ++ ){
            if(subdomain_indices[i] == BACKGROUND_GRID_VALUE) off = 1;
            max_si = std::max( (int)getSubdomainIndex(subdomain_indices[i]), (int)max_si );
        }

        for(unsigned int i = 0 ; i < subdomain_indices.size() ; i++){
            GRID_TYPE si = subdomain_indices[i];
            QColor color(0,0,0);
            if( si > BACKGROUND_GRID_VALUE ){
                color.setHsvF(0.98*(double)(getSubdomainIndex(si)-off)/max_si, 0.8,1.);
                // std::cout << (double)getSubdomainIndex(si)/max_si << " "<< max_si << std::endl;
                if( fileColors.size() > 0 )
                    color = fileColors[getSubdomainIndex(si)];
            }
            colors[i] = color;
        }

#if 0
        int si_size = subdomain_indices.size() ;
        if( counter < 0 ) counter = subdomain_indices.size() - 1;
        if( counter >= subdomain_indices.size() ) counter = 0;
        colors.clear();
        colors.resize(subdomain_indices.size(), QColor(0,0,0) );
        
        std::vector< std::pair<GRID_TYPE, int> > sorted_ids;
        int off = 0;
        for( unsigned int i = 0 ; i < subdomain_indices.size() ; i ++ ){
            if(subdomain_indices[i] == BACKGROUND_GRID_VALUE) off = 1;
            sorted_ids.push_back(std::make_pair(sortedVoxels[i].size(), i));
        }
        
        std::sort( sorted_ids.begin(), sorted_ids.end() );
        
        double step = 1./si_size;

        int s = 0;
        for(unsigned int i = 0 ; i < sorted_ids.size() ; i++){
            
            double h = ((counter+s)%si_size)*step + i%3/3.;

            if( h > 1. ) h -= 1.;
            //h = 1.- h;
            GRID_TYPE si = sorted_ids[i].first;
            QColor color(0,0,0);
            // if( si > BACKGROUND_GRID_VALUE )
            //          color.setHsvF(h/360., 1.,1.);
            color.setHsvF(0.98*h, 0.8,1.);
            
            colors[sorted_ids[i].second] = color;
            //            if(getSubdomainIndex(sorted_ids[i].first) == 48 )
            //                colors[sorted_ids[i].second] = QColor(37,111,253);
            //        if(getSubdomainIndex(sorted_ids[i].first) == 58 )
            //            colors[sorted_ids[i].second] = QColor(202,15,155);
            if(i%3 == 2) s++;
        }
        
        
        //std::cout << sorted_ids.size() << std::endl;
        //colors[1] = QColor(146, 0, 255);
        //std::cout << " count "<< counter << std::endl;
#endif
    }
    
    void VoxelGrid::keepVisiblePart( VoxelGrid & grid ){

        BasicPoint from (cut);
        BasicPoint to ( xdim(), ydim(), zdim() );

        for( int i = 0 ; i < 3 ; i ++ ){
            from[i] = std::max( from[i], (float)0 );
            if( cutDirection[i] < 0 ){
                from[i] = 0;
                to[i] = std::min( cut[i], (float)dim(i) );
            }
        }

        grid.clearAll();

        std::cout << "keepVisiblePart "<< from << " to " << to << std::endl;

        BasicPoint new_dim = to - from;
        grid = VoxelGrid( new_dim[0], new_dim[1], new_dim[2], d(0), d(1), d(2), true, _offset, VOffset );
        grid.clearToDefault();
        for( int i = from[0], ni = 0 ; i < to[0] ; i++, ni++ ){
            for( int j = from[1], nj = 0 ; j <  to[1] ; j++, nj++ ){
                for( int k = from[2], nk = 0 ; k < to[2] ; k++, nk++ ){

                    int id = getGridIndex(i,j,k);
                    GRID_TYPE label = grid_data[id];
                    if(label > BACKGROUND_GRID_VALUE){
                        grid.setValue(ni,nj,nk, label);
                    }
                }
            }
        }

//        grid.buildForDisplay(); TODO check
    }



    void VoxelGrid::removeUnconnectedVoxels( ){

        std::cout << "removeUnconnectedVoxels" << std::endl;

        std::vector<GRID_TYPE> data (grid_data) ;

        for( unsigned int s = 0 ; s < sortedVoxels.size() ;  s++ ){

            std::vector<int> & current_voxels = sortedVoxels[ s ];

            for( unsigned int i = 0 ; i < current_voxels.size(); i++ ){

                int v_id = current_voxels[i];
                Voxel voxel = voxels[v_id];

                int current_grid_index = getGridIndex( voxel );
                GRID_TYPE current_value = value( current_grid_index );

                std::vector<Voxel> neighbors;
                collectSixNeighborhoodVoxels(voxel, neighbors);

                bool keep = false;
                std::map<GRID_TYPE, int> n_ids;

                for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                    int n_id = getGridIndex(neighbors[n]);
                    GRID_TYPE n_grid_value = value(n_id);
                    if( n_grid_value == current_value ){
                        keep = true;
                    }
                    n_ids[ n_grid_value ] ++;
                }

                if( !keep ){
                    int max_v = 0;
                    GRID_TYPE new_value;
                    for( std::map<GRID_TYPE, int>::iterator it = n_ids.begin(); it != n_ids.end(); ++it ){
                        if( it->second > max_v ){
                            max_v = it->second;
                            new_value = it->first;
                        }
                    }
                    data[current_grid_index] = new_value;
                }
            }
        }

        *this = VoxelGrid( data, dim(0), dim(1), dim(2), d(0), d(1), d(2), true, _offset, VOffset );
    }


    void VoxelGrid::fitToData( int offset , VoxelGrid & grid ){

        grid.clear();

        Vmax = Voxel( 0, 0, 0 );
        Vmin = Voxel( _dim[0], _dim[1], _dim[2] );

        for( unsigned int s = 0 ; s < sortedVoxels.size() ;  s++ ){
            std::vector<int> & current_voxels = sortedVoxels[ s ];
            for( unsigned int i = 0 ; i < current_voxels.size(); i++ ){
                int v_id = current_voxels[i];
                Voxel voxel = voxels[v_id];

                Vmax = Voxel::max( Vmax, voxel );
                Vmin = Voxel::min( Vmin, voxel );
            }
        }

        Voxel new_dim = Vmax - Vmin +Voxel( 1, 1, 1 ) ;

        new_dim = Voxel( new_dim.i() + 2*offset, new_dim.j() + 2*offset, new_dim.k() + 2*offset );

        std::cout << "Vmin " << Vmin - Voxel( offset, offset, offset ) << std::endl;

        grid = VoxelGrid( new_dim.i(), new_dim.j(), new_dim.k(), d(0), d(1), d(2), true, _offset, VOffset );
        grid.clearToDefault();

        for( int i = Vmin.i(), ni = offset ; i <= Vmax.i() ; i++, ni++ ){
            for( int j = Vmin.j(), nj = offset ; j <=  Vmax.j() ; j++, nj++ ){
                for( int k = Vmin.k(), nk = offset ; k <= Vmax.k() ; k++, nk++ ){
                    int id = getGridIndex(i,j,k);
                    GRID_TYPE label = grid_data[id];
                    if(label > BACKGROUND_GRID_VALUE){
                        grid.setValue(ni,nj,nk, label);
                    }
                }
            }
        }

//        grid.buildForDisplay();
    }

    void VoxelGrid::addLayer( int id , VoxelGrid & grid ){

        grid.clearAll();

        std::cout << "add layer of " << id  << std::endl;

        std::vector<GRID_TYPE> data (grid_data) ;

        GRID_TYPE new_value = GRID_TYPE(id);
        for(unsigned int i = 0 ; i < _dim[0] ; i++){
            for(unsigned int j = 0 ; j < _dim[1] ; j++){
                for(unsigned int k = 0 ; k < _dim[2] ; k++){
                    int v_id = getGridIndex(i,j,k);
                    GRID_TYPE curren = grid_data[v_id];

                    if( curren <= BACKGROUND_GRID_VALUE ){
                        Voxel voxel = Voxel(i,j,k);

                        Voxel v_born_min, v_born_max;
                        getNeighborhoodBorns( i, j,k, v_born_min, v_born_max, 1);


                        for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                            for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                                for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                    if(vi != voxel.i() || vj != voxel.j() || vk != voxel.k()){
                                        int n_id = getGridIndex(vi,vj,vk);

                                        GRID_TYPE n_value = grid_data[n_id];

                                        if( n_value > BACKGROUND_GRID_VALUE ){
                                            // std::cout << "Separate " << getSubdomainIndex(curren) <<  " and " << getSubdomainIndex(n_value)<< std::endl;
                                            curren = new_value;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    data[ v_id ] = curren;
                }
            }
        }

        grid = VoxelGrid( data, dim(0), dim(1), dim(2), d(0), d(1), d(2), true, _offset, VOffset );

    }

    void VoxelGrid::replaceVisibleLabelsValue( std::map<Subdomain_index, bool> & labels_to_decompose, Subdomain_index to, VoxelGrid & grid ){

        grid.clearAll();

        std::cout << "replaceVisibleLabelsValue to " << to  << std::endl;

        std::vector<GRID_TYPE> data (grid_data) ;

        GRID_TYPE new_value = getGridTypeValue( to );
        for( unsigned int s = 0 ; s < sortedVoxels.size() ;  s++ ){

            std::vector<int> & current_voxels = sortedVoxels[ s ];

            if( labels_to_decompose[s] ){
                std::cout << "Replace label " << s  << " by " << to << std::endl;
                for( unsigned int i = 0 ; i < current_voxels.size(); i++ ){
                    int current_grid_index = getGridIndex( voxels[ current_voxels[i] ] );
                    data[ current_grid_index ] = new_value;
                }
            }
        }

        grid = VoxelGrid( data, dim(0), dim(1), dim(2), d(0), d(1), d(2), true, _offset, VOffset );
    }


    void VoxelGrid::separateInConnectedComponents( std::map<Subdomain_index, bool> & labels_to_decompose, VoxelGrid & grid ){

        grid.clearAll();

        std::cout << "separateInConnectedComponents" << std::endl;

        std::vector<GRID_TYPE> data (grid_data) ;
        std::vector<bool> visited( _size, false );

        GRID_TYPE new_value = BACKGROUND_GRID_VALUE +1;
        for( unsigned int s = 0 ; s < sortedVoxels.size() ;  s++ ){

            std::vector<int> & current_voxels = sortedVoxels[ s ];

            if( labels_to_decompose[s] ){
                std::cout << " Decompose label " << s  << " nb voxels " << current_voxels.size() << std::endl;
                for( unsigned int i = 0 ; i < current_voxels.size(); i++ ){

                    int v_id = current_voxels[i];
                    Voxel voxel = voxels[v_id];

                    int current_grid_index = getGridIndex( voxel );

                    GRID_TYPE current_value = value( current_grid_index );

                    std::vector<Voxel> new_component;

                    if( !visited[ current_grid_index ] ){

                        new_component.clear();
                        std::queue<Voxel> Q;

                        Q.push( voxel );
                        while( std::find( subdomain_indices.begin(), subdomain_indices.end(), new_value ) != subdomain_indices.end() && new_value < CHAR_MAX ){
                            new_value++;
                        }

                        while ( !Q.empty() ){

                            Voxel current_voxel = Q.front();
                            Q.pop();

                            int current_vid = getGridIndex(current_voxel);
                            data[ current_vid ] = new_value;
                            visited[ current_vid ] = true;
                            new_component.push_back( current_voxel );

                            std::vector<Voxel> neighbors;
                            collectSixNeighborhoodVoxels(current_voxel, neighbors);

                            for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                                int n_id = getGridIndex(neighbors[n]);
                                GRID_TYPE n_grid_value = value(n_id);
                                if( n_grid_value == current_value && ! visited[ n_id ]){
                                    Q.push( neighbors[n] );
                                    visited[ n_id ] = true;
                                }
                            }

                        }
                    }


                    subdomain_indices.push_back( new_value );

                }
            }
        }

        grid = VoxelGrid( data, dim(0), dim(1), dim(2), d(0), d(1), d(2), true, _offset, VOffset );
    }

    void VoxelGrid::keepLabelsValues( std::map<Subdomain_index, bool> & label_to_keep, VoxelGrid & grid ){

        grid.clearAll();

        std::map<GRID_TYPE, bool> grid_values_to_keep;

        for( std::map<Subdomain_index, bool>::iterator it = label_to_keep.begin() ; it != label_to_keep.end(); it++ ){
            if( it->second ){
                grid_values_to_keep[ subdomain_indices[ it->first ] ] = true;
            }
        }

        std::vector<GRID_TYPE> data ( _size, BACKGROUND_GRID_VALUE ) ;

        GRID_TYPE default_value = BACKGROUND_GRID_VALUE +1;
        while( std::find( subdomain_indices.begin(), subdomain_indices.end(), default_value ) != subdomain_indices.end() && default_value < CHAR_MAX ){
            default_value++;
        }

        for(int i = 0 ; i < _size ; i ++ ){
            GRID_TYPE value = grid_data[i];
            if( value > BACKGROUND_GRID_VALUE ){
                if( grid_values_to_keep.find( value ) == grid_values_to_keep.end() )
                    value = default_value;

                data[i] = value;
            }
        }

        grid = VoxelGrid( data, dim(0), dim(1), dim(2), d(0), d(1), d(2), true, _offset, VOffset );
    }


    void VoxelGrid::dilateToRemoveVisible( std::map<Subdomain_index, bool> & label_to_remove, VoxelGrid & grid ){
        grid.clearAll();

        std::cout << "dilateToRemoveVisible" << std::endl;

        std::vector<GRID_TYPE> data (grid_data) ;
        std::vector<bool> visited ( grid_data.size(), false );
        for( unsigned int s = 0 ; s < sortedVoxels.size() ;  s++ ){

            std::vector<int> & current_voxels = sortedVoxels[ s ];

            if( label_to_remove[s] ){
                std::cout << " Remove label " << s  << " nb voxels " << current_voxels.size() << std::endl;
                GRID_TYPE current_value = subdomain_indices[s];

                std::queue<Voxel> Q;

                for( unsigned int i = 0 ; i < current_voxels.size(); i++ ){

                    
                    int v_id = current_voxels[i];
                    Voxel voxel = voxels[v_id];

                    std::vector<Voxel> neighbors;
                    collectSixNeighborhoodVoxels(voxel, neighbors);

                    for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                        int n_id = getGridIndex(neighbors[n]);
                        GRID_TYPE n_grid_value = value(n_id);
                        if( n_grid_value != current_value ){
                            Q.push( neighbors[n] );
                        }
                    }
                }

                while( !Q.empty() ){

                    Voxel voxel = Q.front();
                    Q.pop();

                    int current_grid_index = getGridIndex( voxel );
                    visited[ current_grid_index ] = true;
                    GRID_TYPE new_value = value( current_grid_index );

                    std::vector<Voxel> neighbors;
                    collectSixNeighborhoodVoxels(voxel, neighbors);

                    for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                        int n_id = getGridIndex(neighbors[n]);
                        GRID_TYPE n_grid_value = value(n_id);

                        if( !visited[n_id] && n_grid_value == current_value ){
                            data[n_id] = new_value;
                            visited[n_id] = true;
                            Q.push( neighbors[n] );
                        }
                    }
                }
            }
        }

        grid = VoxelGrid( data, dim(0), dim(1), dim(2), d(0), d(1), d(2), true, _offset, VOffset );

    }

    bool VoxelGrid::changeLabelValue( Subdomain_index from, Subdomain_index to ){

        if( from == to ) return false;

        GRID_TYPE si_1 = (GRID_TYPE)from + BACKGROUND_GRID_VALUE;
        GRID_TYPE si_2 = (GRID_TYPE)to + BACKGROUND_GRID_VALUE;

        int si_1_index = -1, si_2_index = -1;
        for( unsigned int i = 0 ; i < subdomain_indices.size(); i++  ){
            if( subdomain_indices[i] == si_1 ){
                si_1_index = i;
            }
            if( subdomain_indices[i] == si_2 ){
                si_2_index = i;
            }
        }

        if( si_1_index >= 0 ){
            if( si_2_index < 0 ){
                si_2_index = sortedVoxels.size();
                sortedVoxels.push_back( std::vector<int>() );
            }
            std::vector<int> & v_indices = sortedVoxels[si_2_index];
            for( unsigned int i = 0 ; i < v_indices.size() ; i++ ){
            }

            return true;
        }

        return false;
    }

    void VoxelGrid::setWhite(){

        for(unsigned int i = 0 ; i < subdomain_indices.size() ; i++){
            colors[i] = Qt::white;
            if( subdomain_indices[i] == BACKGROUND_GRID_VALUE )
                colors[i] = QColor(255,140,0);
        }

    }

    void VoxelGrid::setColors(const std::vector<QColor> &_colors){

        if( colors.size() != _colors.size() ) return;

        colors = _colors;

    }

    void VoxelGrid::getDistanceFiels(std::vector<float> &d_field){
        d_field.clear();
        d_field.resize(_size, 0.);
        std::vector<bool> visited ( _size, false );

        float max = -FLT_MAX;

        for(int i = 0 ; i < _dim[0] ; i++){
            for(int j = 0 ; j <  _dim[1] ; j++){
                for(int k = 0 ; k < _dim[2] ; k++){
                    int g_id = getGridIndex(i,j,k);
                    if( value( g_id ) <= BACKGROUND_GRID_VALUE ){
                        visited[ g_id ] = true;
                    }
                }
            }
        }

        std::priority_queue < std::pair<float, Voxel>, std::deque< std::pair<float, Voxel> > , std::greater< std::pair<float, Voxel> > > DistanceQueue;

        for( unsigned int i = 0 ; i < envelopVoxels.size() ; i ++ ){
            DistanceQueue.push( std::make_pair(10., envelopVoxels[i] ) );
        }

        while( !DistanceQueue.empty()){
            std::pair<float, Voxel> dist_voxel= DistanceQueue.top();
            DistanceQueue.pop();

            float dist = dist_voxel.first;
            Voxel voxel = dist_voxel.second;

            int i = voxel.i();
            int j = voxel.j();
            int k = voxel.k();

            int g_id = getGridIndex( i,j,k );

            if( ! visited[ g_id ] ){

                visited[ g_id ] = true;
                d_field[ g_id ] = dist;

                Voxel v_born_min, v_born_max;
                getNeighborhoodBorns( voxel, v_born_min, v_born_max, 1);

                for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                    for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                        for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                            if(vi != i || vj != j || vk != k){

                                int n_id = getGridIndex( vi, vj, vk );
                                if( ! visited[ n_id ] ){
                                    Voxel n(vi, vj, vk);
                                    float n_dist = dist + ( n - voxel ).norm();
                                    DistanceQueue.push(std::make_pair( n_dist, n ) );
                                }
                            }
                        }
                    }
                }
            }
        }

    }


    bool VoxelGrid::ajustSizeToOffset( ){

        bool sizeChanged = false;

        voffset = Voxel (0,0,0);

        bool update = false;
        for( int i = 0 ; i < 3 ; i++ ){
            int diff = dim(i) - (Vmax[i] + VOffset + 1) ;
            if( diff < 0 )
                update = true;

            diff = Vmin[i] - VOffset;
            if( diff < 0 ){
                voffset[i] += abs(diff);
                update = true;
                sizeChanged = true;
            }
        }

        if(update){
            unsigned int new_dim [] = { _dim[0], _dim[1], _dim[2] };
            unsigned int current_dim [] = { _dim[0], _dim[1], _dim[2] };
            for( int i = 0 ; i < 3 ; i++ ){
                int diff = dim(i) - (Vmax[i] + VOffset + 1) ;
                if( diff < 0 )
                    new_dim[i] -= diff;

                diff = Vmin[i] - VOffset;
                if( diff < 0 ){
                    new_dim[i] += abs(diff);
                }
            }

            if( sortedVoxels.size() > 0 ){

                std::vector<GRID_TYPE> data(voxels.size(), DEFAULT_GRID_VALUE);
                for( unsigned int i = 0 ; i < voxels.size() ; i ++ ){
                    data[i] = value(voxels[i]);
                    voxels[i] += voffset;
                }

                for( unsigned int i = 0 ; i < envelopVoxels.size() ; i ++ ){
                    envelopVoxels[i] += voffset;
                }
                //clearVBOs();

                init(new_dim[0], new_dim[1], new_dim[2], _d[0], _d[1], _d[2]);

                for( unsigned int v = 0 ; v < voxels.size() ; v++ ){
                    setValue(voxels[v], data[v]);
                }

            } else {
                std::vector<GRID_TYPE> current_state(grid_data);

                init(new_dim[0], new_dim[1], new_dim[2], _d[0], _d[1], _d[2]);

                for(int i = 0 ; i < current_dim[0] ; i++){
                    for(int j = 0 ; j <  current_dim[1] ; j++){
                        for(int k = 0 ; k < current_dim[2] ; k++){

                            int current_index = i + j*current_dim[0] + k*current_dim[0]*current_dim[1];
                            int new_index = getGridIndex( i + voffset[0] , j + voffset[1] , k + voffset[2] );
                            setValue( new_index, current_state[current_index]);
                        }
                    }
                }

            }
        }

        return sizeChanged;
        // _offset -= getWorldCoordinate(voffset);
    }



    //void VoxelGrid::reduceSizeToOffset( ){
    //
    //    voffset = Voxel (0,0,0);
    //
    //    Voxel diff = Vmax - Vmin;
    //
    //    std::vector<GRID_TYPE> data(diff[0]*diff[1]*diff[2], 0);
    //       for( unsigned int i = 0 ; i < voxels.size() ; i ++ ){
    //           const Voxel & voxel = voxels[i];
    //           int id = voxel.i() + voxel.j()*diff[0] + voxel.k()*diff[0]*diff[1];
    //            data[i] = value(voxel);
    //            voxels[i] = voffset;
    //        }
    //
    //        for( int i = 0 ; i < 3 ; i++ ){
    //            int diff = dim(i) - (Vmax[i] + VOffset + 1) ;
    //            if( diff < 0 )
    //                _dim[i] -= diff;
    //
    //            diff = Vmin[i] - VOffset;
    //            if( diff < 0 ){
    //                _dim[i] += abs(diff);
    //            }
    //        }
    //
    //        clearVBOs();
    //
    //        init(_dim[0], _dim[1], _dim[2], _d[0], _d[1], _d[2]);
    //
    //        for( unsigned int v = 0 ; v < voxels.size() ; v++ ){
    //            setValue(voxels[v], data[v]);
    //        }
    //
    //    voffset = Voxel (0,0,0);
    //}

    bool VoxelGrid::dilatation( const std::vector<int> & envelopVoxels ){

        GRID_TYPE l_value = BACKGROUND_GRID_VALUE;


        bool sizeChanged = ajustSizeToOffset();

        std::vector<int> next_ids;
        std::vector<int> voxels_ids = envelopVoxels;

        for(unsigned int d = 0 ; d < VOffset; d++){
            for( unsigned int v = 0 ; v < voxels_ids.size() ; v++ ){

                const Voxel & voxel = voxels[voxels_ids[v]];
                std::vector<Voxel> toCheck;

                collectSixNeighborhoodVoxels(voxel, toCheck);

                for(unsigned int l = 0 ; l < toCheck.size() ; l ++){
                    Voxel & vv = toCheck[l];
                    int id_v = getGridIndex(vv);
                    int si_v = grid_data[id_v];

                    if(si_v == BACKGROUND_GRID_VALUE){
                        grid_data[id_v] = l_value;
                        sortedVoxels[0].push_back(voxels.size());
                        next_ids.push_back(voxels.size());
                        voxels.push_back(vv);
                    }
                }
            }

            if( next_ids.size() > 0 ){
                voxels_ids = std::vector<int>(next_ids);
            }
        }

        return sizeChanged;

    }

    //void VoxelGrid::updateAndAddToUsedVoxels(Voxel & voxel){
    //    Voxel & vv = toCheck[l];
    //    int id_v = getGridIndex(vv);
    //    int si_v = grid_data[id_v];
    //
    //    if(si_v < 0){
    //        vv.setLabel(label);
    //        grid_data[id_v] = label;
    //        sortedVoxels[label].push_back(voxels.size());
    //        voxels.push_back(vv);
    //    }
    //
    //}


    void VoxelGrid::computePositions(){

        positions.clear();
        positions.resize(voxels.size());
        for(unsigned int i = 0 ; i < voxels.size() ; i++){
            const Voxel & voxel = voxels[i];
            positions[i] = getWorldCoordinate(voxel);
        }

        computeBoundingBox();
    }

    float VoxelGrid::computeMinPointDist(){

        float min_dist = FLT_MAX;
        for( unsigned int i = 0 ; i < positions.size() ; i ++ ){
            BasicPoint & position = positions[i] ;
            for( unsigned int j = 0 ; j < positions.size() ; j ++ ){
                BasicPoint & neighbor = positions[j];

                min_dist = std::min((neighbor - position).norm() , min_dist);
            }
        }

        return min_dist;
    }

    void VoxelGrid::computeBoundingBox(){

        BBMin = BasicPoint(FLT_MAX, FLT_MAX, FLT_MAX);
        BBMax = BasicPoint(-FLT_MAX, -FLT_MAX, -FLT_MAX);

        for(unsigned int i = 0 ; i < positions.size() ; i++){
            const BasicPoint & point = positions[i];
            for( int v = 0 ; v < 3 ; v++ ){
                BBMin[v] = std::min(BBMin[v], point[v]);
                BBMax[v] = std::max(BBMax[v], point[v]);
            }
        }

    }

    void VoxelGrid::getDeformedVoxelGrid(VoxelGrid & voxelGrid, double resolution){

        std::cout<< "resolution " << resolution << std::endl;
        voxelGrid.clearAll();

        computeBoundingBox();

        BBMax += BasicPoint(0.5,0.5,0.5);
        BBMin -= BasicPoint(0.5,0.5,0.5);

        BasicPoint point = BBMax - BBMin ;

        int n[3];

        double d[3] = {_d[0]*resolution, _d[1]*resolution, _d[2]*resolution};

        for(int i = 0 ; i < 3 ; i++){
            n[i] = int(fabs(point[i])/d[i])+1;
            //        if( n[i] < _dim[i] ){
            //            n[i] = _dim[i];
            //            d[i] = fabs(point[i]/_dim[i]);
            //        }
        }

        // std::cout << n[0] << " " << n[1] << " " << n[2] << std::endl ;

        voxelGrid = VoxelGrid(n[0], n[1], n[2], d[0], d[1], d[2], true, BBMin);

    }

    void VoxelGrid::getDeformedVoxelGrid(VoxelGrid & voxelGrid, const BasicPoint & bbmin, const BasicPoint & bbmax){

        voxelGrid.clearAll();

        BasicPoint point = bbmax - bbmin ;

        int n[3];

        double d[3] = {_d[0],_d[1], _d[2]};

        for(int i = 0 ; i < 3 ; i++){
            n[i] = int(fabs(point[i])/d[i])+1;
        }

        voxelGrid = VoxelGrid(n[0], n[1], n[2], d[0], d[1], d[2], false, bbmin);

    }


    void VoxelGrid::getSegmentation (VoxelGrid & segmentation){

        segmentation = VoxelGrid( dim(0), dim(1), dim(2), dx(), dy(), dz(), false);
        segmentation.clearToZero();
        std::vector<Voxel> & env = segmentation.getEnvelopVoxels();
        for(int i = 0 ; i < envelopVoxels.size() ; i++)
            env.push_back(envelopVoxels[i]);

        for(int i = 0 ; i < _dim[0] ; i++){
            for(int j = 0 ; j <  _dim[1] ; j++){
                for(int k = 0 ; k < _dim[2] ; k++){

                    int id = getGridIndex(i,j,k);
                    GRID_TYPE label = BACKGROUND_GRID_VALUE;
                    if(grid_data[id] > BACKGROUND_GRID_VALUE)
                        label += 1;
                    segmentation.setValue(i,j,k,label);
                }
            }
        }
    }

    bool VoxelGrid::isPairToSeparate( Subdomain_index l1, Subdomain_index l2 ){

        if( l1 == 4 && l2 == 9)
            return true;
        if( l1 == 4 && l2 == 10)
            return true;
        if( l1 == 4 && l2 == 11)
            return true;
        if( l1 == 5 && l2 == 11)
            return true;
        if( l1 == 5 && l2 == 10)
            return true;
        if( l1 == 7 && l2 == 9)
            return true;
        if( l1 == 7 && l2 == 10)
            return true;
        if( l1 == 7 && l2 == 14)
            return true;
        if( l1 == 8 && l2 == 10)
            return true;
        if( l1 == 8 && l2 == 14)
            return true;
        return false;
    }

    void VoxelGrid::automaticSeparation(){

        std::vector<GRID_TYPE> ids_to_sep(7);
        ids_to_sep[0] = getGridTypeValue(4); ids_to_sep[1] = getGridTypeValue(5); ids_to_sep[2] = getGridTypeValue(7);
        ids_to_sep[3] = getGridTypeValue(8); ids_to_sep[4] = getGridTypeValue(10); ids_to_sep[5] = getGridTypeValue(11); ids_to_sep[6] = getGridTypeValue(14);

        std::vector<Voxel> borderV ;
        std::vector<GRID_TYPE> current_state = std::vector<GRID_TYPE> (grid_data);
        for(unsigned int i = 0 ; i < _dim[0] ; i++){
            for(unsigned int j = 0 ; j < _dim[1] ; j++){
                for(unsigned int k = 0 ; k < _dim[2] ; k++){
                    int v_id = getGridIndex(i,j,k);
                    GRID_TYPE curren = current_state[v_id];

                    if( curren > BACKGROUND_GRID_VALUE && std::find( ids_to_sep.begin(), ids_to_sep.end(), curren ) != ids_to_sep.end() ){
                        Voxel voxel = Voxel(i,j,k);

                        Voxel v_born_min, v_born_max;
                        getNeighborhoodBorns( i, j,k, v_born_min, v_born_max, 1);


                        for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                            for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                                for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                    if(vi != voxel.i() || vj != voxel.j() || vk != voxel.k()){
                                        int n_id = getGridIndex(vi,vj,vk);

                                        GRID_TYPE n_value = current_state[n_id];

                                        if( n_value > BACKGROUND_GRID_VALUE && curren != n_value && isPairToSeparate( getSubdomainIndex(curren) , getSubdomainIndex(n_value) ) ){
                                            // std::cout << "Separate " << getSubdomainIndex(curren) <<  " and " << getSubdomainIndex(n_value)<< std::endl;
                                            curren = DEFAULT_GRID_VALUE;
                                            borderV.push_back(Voxel(vi,vj,vk));
                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        setValue(v_id, curren);

                    }
                }
            }
        }




        for( int v = 1 ; v < abs(15) ; v++ ){
            current_state.clear();
            std::vector<GRID_TYPE> current_state = std::vector<GRID_TYPE> (grid_data);
            std::vector<Voxel> next;
            for( unsigned int b = 0 ; b < borderV.size() ; b++ ){

                int i = borderV[b].i();
                int j = borderV[b].j();
                int k = borderV[b].k();
                int v_id = getGridIndex(i,j,k);

                Voxel voxel = Voxel(i,j,k);

                setValue(v_id, DEFAULT_GRID_VALUE);
                GRID_TYPE curren = current_state[v_id];
                Voxel v_born_min, v_born_max;
                getNeighborhoodBorns( i, j,k, v_born_min, v_born_max, 1);
                if(curren > BACKGROUND_GRID_VALUE)
                    for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                        for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                            for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                if(vi != voxel.i() || vj != voxel.j() || vk != voxel.k()){
                                    setValue(vi, vj,vk, BACKGROUND_GRID_VALUE);
                                    next.push_back(Voxel(vi,vj,vk));
                                }
                            }
                        }
                    }
            }
            borderV = next;
        }

        if( userOffsets.size() == 0 )
            userOffsets.resize( envelopVoxels.size(), 1 );
        for( unsigned int b = 0 ; b < borderV.size() ; b++ ){
            envelopVoxels.push_back(borderV[b]);
            userOffsets.push_back(-1);
        }

    }

    void VoxelGrid::computeEnvelopVoxels(){

        envelopVoxels.clear();

        for(int i = 0 ; i < _dim[0] ; i++){
            for(int j = 0 ; j <  _dim[1] ; j++){
                for(int k = 0 ; k < _dim[2] ; k++){

                    int id = getGridIndex(i,j,k);
                    GRID_TYPE label = grid_data[id];

                    if(label > BACKGROUND_GRID_VALUE){
                        std::vector<Voxel> neighbors;
                        collectSixNeighborhoodVoxels(i, j, k, neighbors);

                        for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                            const Voxel & voxel = neighbors[n];

                            if( value( voxel ) <= BACKGROUND_GRID_VALUE ){
                                envelopVoxels.push_back(Voxel(i,j,k));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    void VoxelGrid::dilatation( int voffsetNb , unsigned int dx, unsigned int dy, unsigned int dz, const Voxel & v_offset, bool s_changed){


        if(voffsetNb == VOffset) return ;



#if 0
        int it_nb = voffsetNb - VOffset;
        bool dil = true;
        if(it_nb < 0 ){
            dil = false;
            it_nb *= -1;
        }
        std::vector<Voxel> borderV ;
        std::vector<GRID_TYPE> current_state = std::vector<GRID_TYPE> (grid_data);

        if(!dil){
            for(unsigned int i = 0 ; i < _dim[0] ; i++){
                for(unsigned int j = 0 ; j < _dim[1] ; j++){
                    for(unsigned int k = 0 ; k < _dim[2] ; k++){
                        int v_id = getGridIndex(i,j,k);
                        GRID_TYPE curren = current_state[v_id];
                        if(curren > BACKGROUND_GRID_VALUE){

                            Voxel voxel = Voxel(i,j,k);

                            Voxel v_born_min, v_born_max;
                            getNeighborhoodBorns( i, j,k, v_born_min, v_born_max, 1);


                            for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                                for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                                    for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                        if(vi != voxel.i() || vj != voxel.j() || vk != voxel.k()){

                                            int n_id = getGridIndex(vi,vj,vk);

                                            GRID_TYPE n_value = current_state[n_id];
                                            if( n_value <= BACKGROUND_GRID_VALUE ){

                                                curren = BACKGROUND_GRID_VALUE;
                                                borderV.push_back(Voxel(vi,vj,vk));
                                                break;
                                            }
                                        }
                                    }
                                }
                            }


                            setValue(v_id, curren);
                        }
                    }
                }
            }

            for( int v = 1 ; v < it_nb ; v++ ){
                current_state.clear();
                std::vector<GRID_TYPE> current_state = std::vector<GRID_TYPE> (grid_data);
                std::vector<Voxel> next;
                for( unsigned int b = 0 ; b < borderV.size() ; b++ ){

                    int i = borderV[b].i();
                    int j = borderV[b].j();
                    int k = borderV[b].k();
                    int v_id = getGridIndex(i,j,k);

                    Voxel voxel = Voxel(i,j,k);

                    setValue(v_id, DEFAULT_GRID_VALUE);
                    GRID_TYPE curren = current_state[v_id];
                    Voxel v_born_min, v_born_max;
                    getNeighborhoodBorns( i, j,k, v_born_min, v_born_max, 1);
                    if(curren > BACKGROUND_GRID_VALUE)
                        for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                            for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                                for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                    if(vi != voxel.i() || vj != voxel.j() || vk != voxel.k()){
                                        setValue(vi, vj,vk, BACKGROUND_GRID_VALUE);
                                        next.push_back(Voxel(vi,vj,vk));
                                    }
                                }
                            }
                        }
                }
                borderV = next;
            }
        } else {

            for(unsigned int i = 0 ; i < _dim[0] ; i++){
                for(unsigned int j = 0 ; j < _dim[1] ; j++){
                    for(unsigned int k = 0 ; k < _dim[2] ; k++){
                        int v_id = getGridIndex(i,j,k);
                        GRID_TYPE curren = current_state[v_id];
                        if(curren > BACKGROUND_GRID_VALUE){

                            Voxel voxel = Voxel(i,j,k);

                            Voxel v_born_min, v_born_max;
                            getNeighborhoodBorns( i, j,k, v_born_min, v_born_max, 1);


                            for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                                for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                                    for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                        if(vi != voxel.i() || vj != voxel.j() || vk != voxel.k()){

                                            int n_id = getGridIndex(vi,vj,vk);

                                            GRID_TYPE n_value = current_state[n_id];
                                            GRID_TYPE up_value = value(n_id);
                                            if( n_value == BACKGROUND_GRID_VALUE && up_value == n_value){
                                                setValue(n_id, curren);
                                                borderV.push_back(Voxel(vi,vj,vk));
                                                //break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for( int v = 1 ; v < it_nb ; v++ ){
                current_state.clear();
                std::vector<GRID_TYPE> current_state = std::vector<GRID_TYPE> (grid_data);
                std::vector<Voxel> next;
                for( unsigned int b = 0 ; b < borderV.size() ; b++ ){

                    int i = borderV[b].i();
                    int j = borderV[b].j();
                    int k = borderV[b].k();
                    int v_id = getGridIndex(i,j,k);

                    Voxel voxel = Voxel(i,j,k);

                    GRID_TYPE curren = current_state[v_id];
                    Voxel v_born_min, v_born_max;
                    getNeighborhoodBorns( i, j,k, v_born_min, v_born_max, 1);

                    for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                        for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                            for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                if(vi != voxel.i() || vj != voxel.j() || vk != voxel.k()){
                                    int n_id = getGridIndex(vi,vj,vk);

                                    GRID_TYPE n_value = current_state[n_id];
                                    GRID_TYPE up_value = value(n_id);
                                    if( n_value == BACKGROUND_GRID_VALUE && up_value == n_value){
                                        setValue(n_id, curren);
                                        next.push_back(Voxel(vi,vj,vk));
                                        //break;
                                    }
                                }
                            }
                        }

                    }
                    borderV = next;
                }

            }

        }
#else 
        int it_nb = voffsetNb - VOffset;

        VOffset = voffsetNb;

        if(s_changed){

            for( unsigned int i = 0 ; i < envelopVoxels.size() ; i ++ ){
                envelopVoxels[i] += v_offset;
            }

            std::vector<GRID_TYPE> current_state(grid_data);

            int current_dim[3] = {dim(0), dim(1), dim(2)};
            init(dx, dy, dz, _d[0], _d[1], _d[2]);

            for(int i = 0 ; i < current_dim[0] ; i++){
                for(int j = 0 ; j <  current_dim[1] ; j++){
                    for(int k = 0 ; k < current_dim[2] ; k++){

                        int current_index = i + j*current_dim[0] + k*current_dim[0]*current_dim[1];
                        int new_index = getGridIndex( i + v_offset[0] , j + v_offset[1] , k + v_offset[2] );
                        setValue( new_index, current_state[current_index]);
                    }
                }
            }

            std::cout << "Size changed " << std::endl;
        }

        std::cout << "Dilatation " << it_nb << std::endl;
        if( userOffsets.size() > 0 ){

            std::vector<bool> marked(_size, false);
            for(unsigned int i = 0 ; i < envelopVoxels.size() ; i++ ){
                marked[getGridIndex(envelopVoxels[i])] = true;
                userOffsets[i] *= fabs(it_nb + 1);
            }

            VOffset = voffsetNb;

            bool to_do = true;

            std::vector<Voxel> currentEnvelop;
            std::vector<int> currentOffsets;

            while(to_do){
                to_do = false;
                for(unsigned int i = 0 ; i < envelopVoxels.size() ; i++ ){

                    Voxel & voxel = envelopVoxels[i];
                    GRID_TYPE label = value(voxel);

                    if( userOffsets[i] > 1 ){
                        std::cout << "Dilatation " << userOffsets[i] << std::endl;
                        std::vector<Voxel> neighbors;
                        collectSixNeighborhoodVoxels(voxel.i(), voxel.j(), voxel.k(), neighbors);

                        for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                            const Voxel & n_voxel = neighbors[n];

                            if(value( n_voxel ) <= BACKGROUND_GRID_VALUE ){
                                setValue(n_voxel, label);

                                if( userOffsets[i]-1 > 1 ){
                                    currentEnvelop.push_back(n_voxel);
                                    currentOffsets.push_back(userOffsets[i]-1);
                                    to_do = true;
                                }
                            }
                        }

                    } else  if ( userOffsets[i] < -1 ) {
                        std::cout << "Erosion " << userOffsets[i] << std::endl;
                        setValue(voxel, BACKGROUND_GRID_VALUE);

                        if( userOffsets[i] + 1 < -1 ){

                            std::vector<Voxel> neighbors;
                            collectSixNeighborhoodVoxels(voxel.i(), voxel.j(), voxel.k(), neighbors);

                            for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                                const Voxel & n_voxel = neighbors[n];

                                int n_id = getGridIndex(n_voxel);
                                if(value( n_voxel ) > BACKGROUND_GRID_VALUE && !marked[n_id] ){
                                    marked[n_id] = true;
                                    currentEnvelop.push_back(n_voxel);
                                    currentOffsets.push_back(userOffsets[i]+1);
                                    to_do = true;

                                }
                            }
                        }
                    }


                }
                envelopVoxels = currentEnvelop;
                userOffsets = currentOffsets;

                currentEnvelop.clear();
                userOffsets.clear();
            }
        }  else {


            VOffset = voffsetNb;

            bool offset = false;
            if(it_nb > 0) offset = true;

            if(! offset ) it_nb *= -1;

            if( envelopVoxels.size() == 0  ){
                for(unsigned int i = 0 ; i < _dim[0] ; i++){
                    for(unsigned int j = 0 ; j < _dim[1] ; j++){
                        for(unsigned int k = 0 ; k < _dim[2] ; k++){

                            int index = getGridIndex(i,j,k);

                            GRID_TYPE v_value = value(index);

                            if(v_value> BACKGROUND_GRID_VALUE){

                                std::vector<Voxel> neighbors;
                                collectSixNeighborhoodVoxels(i, j, k, neighbors);

                                for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                                    const Voxel & n_voxel = neighbors[n];

                                    if( grid_data [getGridIndex( n_voxel )] <= BACKGROUND_GRID_VALUE || i ==0 || j ==0 || k ==0 || i == _dim[0]-1 || j == _dim[0]-1 || k == _dim[0]-1   ){
                                        envelopVoxels.push_back(Voxel(i,j,k));
                                        break;
                                    }
                                }

                            }
                        }
                    }
                }

            }

            for(unsigned int it = 0 ; it < it_nb ; it++ ){


                std::vector<Voxel> currentEnvelop;

                if ( offset ){

                    for(unsigned int i = 0 ; i < envelopVoxels.size() ; i++ ){
                        Voxel & voxel = envelopVoxels[i];
                        GRID_TYPE label = value ( voxel );


                        std::vector<Voxel> neighbors;
                        collectSixNeighborhoodVoxels(voxel.i(), voxel.j(), voxel.k(), neighbors);

                        for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                            const Voxel & n_voxel = neighbors[n];

                            if( value( n_voxel ) == BACKGROUND_GRID_VALUE ){
                                setValue(n_voxel, label);
                                currentEnvelop.push_back(n_voxel);
                            }
                        }
                    }

                } else{
                    std::vector<bool> marked(_size, false);
                    for(unsigned int i = 0 ; i < envelopVoxels.size() ; i++ ){
                        Voxel & voxel = envelopVoxels[i];

                        setValue(voxel, BACKGROUND_GRID_VALUE);

                        std::vector<Voxel> neighbors;
                        collectSixNeighborhoodVoxels(voxel.i(), voxel.j(), voxel.k(), neighbors);

                        for( unsigned int n = 0 ; n < neighbors.size() ; n++ ){
                            const Voxel & n_voxel = neighbors[n];
                            int n_id = getGridIndex(n_voxel);
                            if( value( n_voxel ) >= BACKGROUND_GRID_VALUE && !marked[n_id]){
                                marked[n_id] = true;
                                currentEnvelop.push_back(n_voxel);
                            }
                        }
                    }
                }

                envelopVoxels = currentEnvelop;
                currentEnvelop.clear();
            }

        }
#endif

    }


    void VoxelGrid::rasterize(VoxelGrid & deformed_grid, double resolution){

        std::map<unsigned int, unsigned int> VMap;
        projectInRegularGrid(deformed_grid, VMap, resolution);

        VMap.clear();

        bool again = true;

        while(again){

            again = false;
#if 1
            const std::vector<GRID_TYPE> currentState(deformed_grid.data());


            for(int i = 0 ; i < deformed_grid.xdim() ; i++){
                for(int j = 0 ; j < deformed_grid.ydim() ; j++){
                    for(int k = 0 ; k < deformed_grid.zdim() ; k++){

                        int id = deformed_grid.getGridIndex(i,j,k);
                        GRID_TYPE label = currentState[id];

                        if(label == BACKGROUND_GRID_VALUE){
                            deformed_grid.dilatate(i, j, k);
                        }
                    }
                }
            }

            for(int i = 0 ; i < deformed_grid.xdim() ; i++){
                for(int j = 0 ; j < deformed_grid.ydim() ; j++){
                    for(int k = 0 ; k < deformed_grid.zdim() ; k++){

                        int id = deformed_grid.getGridIndex(i,j,k);
                        GRID_TYPE label = currentState[id];

                        if(label > BACKGROUND_GRID_VALUE){
                            if(deformed_grid.dilatate(i,j,k))
                                again = true;
                        }
                    }
                }
            }
#else
            again = false;

            const std::vector<GRID_TYPE> currentState = std::vector<GRID_TYPE> (deformed_grid.data());
            std::vector<Voxel> nextStep;


            for(int i = 0 ; i < deformed_grid.xdim() ; i++){
                for(int j = 0 ; j < deformed_grid.ydim() ; j++){
                    for(int k = 0 ; k < deformed_grid.zdim() ; k++){

                        int id = deformed_grid.getGridIndex(i,j,k);
                        int label = currentState[id];

                        if(label == 0){
                            deformed_grid.dilatate(i, j, k);
                        } else if(label > 0){
                            nextStep.push_back(Voxel(i,j,k));
                        }
                    }
                }
            }

            for(int v = 0 ; v < nextStep.size() ; v++){

                int id = deformed_grid.getGridIndex(nextStep[v]);
                int label = currentState[id];

                if(label > 0){
                    if(deformed_grid.dilatate(nextStep[v]))
                        again = true;
                }

            }




            //        if(deformed_grid.dilatate())
            //            again = true;

#endif
        }


        deformed_grid.sort();

        //    deformed_grid = VoxelGrid(tmp, deformed_grid.xdim(), deformed_grid.ydim(), deformed_grid.zdim(),
        //                              deformed_grid.dx(), deformed_grid.dy(), deformed_grid.dz(), deformed_grid.offset());




    }

    void VoxelGrid::projectInRegularGrid(VoxelGrid & deformed_grid, std::map<unsigned int, unsigned int> & VMap, double resolution, bool sort){

        getDeformedVoxelGrid(deformed_grid, resolution);

        for(unsigned int i = 0 ; i < deformed_grid.size() ; i++){
            deformed_grid.setValue( i, DEFAULT_GRID_VALUE );
        }


        for(unsigned int v = 0; v < voxels.size(); v++){
            Voxel projected_voxel = deformed_grid.getGridCoordinate(positions[v]);
            unsigned int id = deformed_grid.getGridIndex(projected_voxel);
            deformed_grid.setValue( id, value(voxels[v]) );
            VMap[id] = v;
        }

        if(sort)
            deformed_grid.sort();
    }

    void VoxelGrid::rasterizeClosest(VoxelGrid & deformed_grid, unsigned int neighborhoodSize, double resolution){

        std::map<unsigned int, unsigned int> VMap;
        projectInRegularGrid(deformed_grid, VMap, resolution);

        std::vector<GRID_TYPE> current_state(deformed_grid.data());
        for(unsigned int i = 0 ; i < deformed_grid.xdim() ; i++ ){
            for(unsigned int j = 0 ; j < deformed_grid.ydim() ; j++ ){
                for(unsigned int k = 0 ; k < deformed_grid.zdim() ; k++ ){
                    unsigned int id = deformed_grid.getGridIndex( i, j, k );

                    if(deformed_grid.value(id) == DEFAULT_GRID_VALUE){
                        Voxel voxel(i,j,k);

                        Voxel v_born_min, v_born_max;
                        deformed_grid.getNeighborhoodBorns( i, j, k, v_born_min, v_born_max, neighborhoodSize);

                        float min_dist = FLT_MAX;
                        GRID_TYPE min_sub_id = DEFAULT_GRID_VALUE;
                        for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                            for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                                for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                    if(vi != i || vj != j || vk != k){

                                        int neig_id = deformed_grid.getGridIndex(vi, vj, vk);
                                        GRID_TYPE si_v = current_state[neig_id];

                                        if(si_v >= BACKGROUND_GRID_VALUE ){
                                            float dist = ( Voxel(vi, vj, vk) - voxel ).norm();
                                            if(dist < min_dist){
                                                min_dist = dist;
                                                min_sub_id = si_v;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if( min_sub_id >=  BACKGROUND_GRID_VALUE)
                            deformed_grid.setValue(i,j,k, min_sub_id);
                    }
                }
            }
        }


        deformed_grid.sort();
    }

    //void VoxelGrid::detectOutliers(){
    //
    //}

    void VoxelGrid::rasterizeUsingVectorFields(VoxelGrid & deformed_grid, unsigned int neighborhoodSize, double resolution){

#if 0
        std::map<unsigned int, unsigned int> VMap;
        projectInRegularGrid(deformed_grid, VMap);

        std::vector<GRID_TYPE> current_state(deformed_grid.data());
        for(unsigned int i = 0 ; i < deformed_grid.xdim() ; i++ ){
            for(unsigned int j = 0 ; j < deformed_grid.ydim() ; j++ ){
                for(unsigned int k = 0 ; k < deformed_grid.zdim() ; k++ ){
                    unsigned int id = deformed_grid.getGridIndex( i, j, k );

                    if(deformed_grid.value(id) < 0){
                        Voxel voxel(i,j,k);

                        Voxel v_born_min, v_born_max;
                        deformed_grid.getNeighborhoodBorns( i, j, k, v_born_min, v_born_max, neighborhoodSize);

                        float min_dist = FLT_MAX;
                        int min_sub_id = -1;

                        float sumW = 0.;
                        int count = 0;

                        for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                            for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                                for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                    if(vi != i || vj != j || vk != k){

                                        int neig_id = deformed_grid.getGridIndex(vi, vj, vk);
                                        int si_v = current_state[neig_id];

                                        if(si_v >= 0 ){
                                            float dist = ( Voxel(vi, vj, vk) - voxel ).norm();
                                            if(dist < min_dist){
                                                min_dist = dist;
                                                min_sub_id = si_v;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if( min_sub_id >=  0)
                            deformed_grid.setValue(i,j,k, min_sub_id);
                    }
                }
            }
        }


        deformed_grid.sort();
#else

        QTime t;
        t.start();

        std::map<unsigned int, unsigned int> VMap;
        projectInRegularGrid(deformed_grid, VMap, resolution);

        std::vector<GRID_TYPE> current_state(deformed_grid.data());
        for(unsigned int i = 0 ; i < deformed_grid.xdim() ; i++ ){
            for(unsigned int j = 0 ; j < deformed_grid.ydim() ; j++ ){
                for(unsigned int k = 0 ; k < deformed_grid.zdim() ; k++ ){
                    unsigned int id = deformed_grid.getGridIndex( i, j, k );

                    if(deformed_grid.value(id) == DEFAULT_GRID_VALUE){
                        Voxel voxel(i,j,k);
                        Voxel v_born_min, v_born_max;
                        deformed_grid.getNeighborhoodBorns( i, j, k, v_born_min, v_born_max, neighborhoodSize );

                        BasicPoint t( 0, 0, 0 );
                        float sumW = 0.;
                        int count = 0;

                        for(unsigned int vi = v_born_min.i() ; vi <= v_born_max.i() ; vi++){
                            for(unsigned int vj = v_born_min.j() ; vj <= v_born_max.j() ; vj++){
                                for(unsigned int vk = v_born_min.k() ; vk <= v_born_max.k() ; vk++){
                                    if(vi != i || vj != j || vk != k){

                                        int neig_id = deformed_grid.getGridIndex(vi, vj, vk);
                                        GRID_TYPE si_v = current_state[neig_id];

                                        if(si_v >= BACKGROUND_GRID_VALUE ){

                                            Voxel projectedNeigh(vi, vj, vk);
                                            Voxel & neigh = voxels[VMap[deformed_grid.getGridIndex(projectedNeigh)]];

                                            float sqrdist = ( Voxel(vi, vj, vk) - voxel ).sqrnorm();

                                            Voxel diff = projectedNeigh- neigh;
                                            float w = 1./sqrdist;
                                            t += w*BasicPoint( diff.i(), diff.j(), diff.k() );

                                            sumW += w;

                                            count ++;
                                        }
                                    }
                                }
                            }
                        }

                        if(count > 0){

                            t = t/sumW;

                            Voxel res =  voxel - Voxel(round(t[0]), round(t[1]), round(t[2]));
                            //  std::cout << res << std::endl;
                            if((res.i() >= 0 && res.j() >= 0 && res.k() >= 0 ) && isInGrid(res))
                                deformed_grid.setValue( voxel, value(res));
                        }

                    }
                }
            }
        }

        std::cout << "Time elapsed for rasterization: " << t.elapsed() << " ms" << std::endl;
        current_state.clear();
        deformed_grid.sort();
#endif
    }

    int VoxelGrid::getVoxelNumber( const Subdomain_index si ){
        std::vector< int >::iterator it = std::find(vnumbers.begin(), vnumbers.end(), getGridTypeValue(si));
        int result = 0;
        if( it != vnumbers.end())
            result = *it;
        return result;
    }

    GRID_TYPE VoxelGrid::getGridTypeValue(Subdomain_index si){
        return si + BACKGROUND_GRID_VALUE;
    }

    Subdomain_index VoxelGrid::getSubdomainIndex(GRID_TYPE g_value){
        return g_value - BACKGROUND_GRID_VALUE;
    }

    void VoxelGrid::getSubdomainIndices(std::vector<Subdomain_index> & _subdomain_indices){
        _subdomain_indices.clear();
        _subdomain_indices.resize(subdomain_indices.size());

        for( unsigned int i = 0 ; i < subdomain_indices.size(); i++)
            _subdomain_indices[i] = getSubdomainIndex(subdomain_indices[i]);
    }

    void VoxelGrid::compareVolumeLossWithGrid( VoxelGrid & grid , std::vector< float > & volume_diff ){

        volume_diff.clear();
        std::vector<int> def_v_nbs;
        std::vector<int> original_v_nbs;
        float min_value = FLT_MAX;
        float max_value = -FLT_MAX;

        int vsum =0;

        float originalVoxelVolume = _d[0]*_d[1]*_d[2];
        float defVoxelVolume = grid.d(0)*grid.d(1)*grid.d(2);
        int max_nb =0;
        int min_nb =0;


        for(unsigned int i = 0 ; i < subdomain_indices.size() ; i++ ){
            GRID_TYPE si = subdomain_indices[i];

            if( si > BACKGROUND_GRID_VALUE ){

                int vNb = getVoxelNumber(si);
                int def_v_nb = grid.getVoxelNumber(si);
                float diff = (def_v_nb*defVoxelVolume - vNb*originalVoxelVolume);

                volume_diff.push_back(diff);

                original_v_nbs.push_back(vNb);
                def_v_nbs.push_back(def_v_nb);

                vsum += vNb;

                if( min_value > diff ){
                    min_nb = vNb;
                    min_value = std::min(diff, min_value);
                }
                if( max_value < diff ){
                    max_nb = vNb;
                    max_value = std::max(diff, max_value);
                }

            }
        }

        float diffPercent = 0.;
        float totalVolumeChange = 0.;
        float originalTotalVolume = vsum*originalVoxelVolume;
        for(unsigned int i = 0 ; i < volume_diff.size() ; i++ ){
            totalVolumeChange += volume_diff[i];
            diffPercent += fabs(volume_diff[i]);
            volume_diff[i] = 100*volume_diff[i]/original_v_nbs[i];
        }

        // volume_diff.push_back(originalTotalVolume);

        std::cout << "Grid of " << vsum << " voxels" << std::endl ;
        std::cout << "Relative max " << 100*max_value/float(max_nb*originalVoxelVolume) << " representing " << 100*float(max_nb)/vsum << " percent of the volume " << std::endl;
        std::cout << "Relative min " << 100*min_value/float(min_nb*originalVoxelVolume) << " representing " <<  100*float(min_nb)/vsum << " percent of the volume " <<  std::endl;
        std::cout << "Volume change " << diffPercent <<  " % "<< std::endl ;
        std::cout << "Volume loss " << totalVolumeChange<< std::endl;
        std::cout <<  "Min " << 100.*min_value/float(originalTotalVolume)<<  ", Max " << 100.*max_value/float(originalTotalVolume) << ",  average " << diffPercent /volume_diff.size() << std::endl;

    }

    void VoxelGrid::computeDistortion(const std::vector<BasicPoint> &defPositions, std::vector<float> &distortions){
        distortions.clear();
        distortions.resize(size(), 0.);


        std::vector<int> VMap(size(), -1);

        for( unsigned int i = 0 ; i < voxels.size() ; i++ ){
            VMap[getGridIndex(voxels[i])] = i;
        }

        for( unsigned int i = 0 ; i < voxels.size() ; i++ ){
            Voxel & voxel = voxels[i];
            int grid_id = getGridIndex(voxel);
            if( value(grid_id) > BACKGROUND_GRID_VALUE ){
                const BasicPoint & defPosition = defPositions[i];

                std::vector<Voxel> neighbors;
                collectSixNeighborhoodVoxels(voxel, neighbors);

                for( unsigned int j = 0 ; j < neighbors.size() ; j++ ){

                    Voxel & neighbor = neighbors[j];
                    unsigned int id = getGridIndex(neighbor);
                    int nid = VMap[id];

                    if( value(id) > BACKGROUND_GRID_VALUE && nid >= 0 ){
                        const BasicPoint & neighborDefPosition = defPositions[nid];

                        BasicPoint defVec = neighborDefPosition - defPosition;
                        defVec.normalize();

                        for( unsigned int k = j+1 ; k < neighbors.size() ; k++ ){

                            Voxel & orthoNeighbor = neighbors[k];
                            unsigned int orthoId = getGridIndex(orthoNeighbor);

                            if( (orthoNeighbor-voxel + neighbor-voxel) != Voxel(0,0,0)){

                                int north = VMap[orthoId];
                                if( value(orthoId) > BACKGROUND_GRID_VALUE  && north >= 0){

                                    BasicPoint orthoDefVec = defPositions[north] - defPosition;
                                    orthoDefVec.normalize();

                                    float angle = acos(dot(defVec, orthoDefVec))*180./M_PI;

                                    distortions[grid_id] = angle;
                                    float dotP = dot(defVec, orthoDefVec);
                                    if( std::isnan(dotP) || std::isnan(angle) ){
                                        std::cout << "dot : " << dotP << std::endl ;
                                        std::cout << "angle : " << acos(dot(defVec, orthoDefVec)) << std::endl ;
                                    }

                                }
                            }
                        }

                    }
                }
            }
        }
        VMap.clear();
    }

    void VoxelGrid::computeDistortion(const std::vector<BasicPoint> & defPositions,  std::vector<int> & values, float & min_value, float & max_value){

        std::vector<float> distortions;
        computeDistortion(defPositions, distortions);

        min_value = 0.;
        max_value = 180;
        int v_nb = 200;

        values.clear();
        values.resize(v_nb, 0);

        float max_angle = -FLT_MAX;
        float min_angle = FLT_MAX;
        float max_concentration = 0;
        int max_id = 0;

        for( unsigned int d = 0 ; d < distortions.size() ; d++ ){
            int j = int(v_nb*(std::max(std::min(distortions[d], max_value),min_value)-min_value)/(max_value-min_value));
            if(j < 0 || j >= values.size())
                std::cout << "distortions " << d << " : " << distortions[d] << " j : " << j << std::endl;
            else {
                if(value(d)>BACKGROUND_GRID_VALUE){
                    values[j] ++;
                    if( distortions[d] > max_angle )
                        max_angle = distortions[d];
                    if( distortions[d] < min_angle )
                        min_angle = distortions[d];

                    if(max_id < values[j]){
                        max_concentration = distortions[d];
                    }
                }
            }

        }


        std::cout << "min : " << min_angle << ", max : " << max_angle << " - " << max_concentration << std::endl;

    }

    void VoxelGrid::computeElongation(const std::vector<BasicPoint> & defPositions, std::vector<int> & values, float & min_value, float & max_value){

        std::vector<float> elongations ;

        std::vector<int> VMap(size(), -1);

        for( unsigned int i = 0 ; i < voxels.size() ; i++ ){
            VMap[getGridIndex(voxels[i])] = i;
        }

        min_value = FLT_MAX;
        max_value = -FLT_MAX;

        for( unsigned int i = 0 ; i < voxels.size() ; i++ ){
            Voxel & voxel = voxels[i];

            if( value(voxel) > BACKGROUND_GRID_VALUE ){
                const BasicPoint & defPosition = defPositions[i];
                BasicPoint position = getWorldCoordinate(voxel);

                std::vector<Voxel> neighbors;
                collectSixNeighborhoodVoxels(voxel, neighbors);

                for( unsigned int j = 0 ; j < neighbors.size() ; j++ ){

                    Voxel & neighbor = neighbors[j];
                    unsigned int id = getGridIndex(neighbor);
                    int nid = VMap[id];

                    if( value(id) > BACKGROUND_GRID_VALUE && nid>=0 ){
                        const BasicPoint & neighborDefPosition = defPositions[nid];
                        BasicPoint neighborPosition = getWorldCoordinate(neighbor);

                        BasicPoint defVec = neighborDefPosition - defPosition;
                        BasicPoint vec = neighborPosition - position;
                        float elongation = defVec.norm()/vec.norm();

                        if (elongation > 3. ) elongation = 1;
                        elongations.push_back(elongation);

                        max_value = std::max(max_value, elongation);
                    }
                }
            }
        }

        VMap.clear();

        int v_nb = 200;

        values.clear();
        values.resize(v_nb, 0);

        min_value = 0.;
        max_value = std::max((float)2., max_value);

        std::cout << "min : " << min_value << ", max : " << max_value << std::endl;

        for( unsigned int i = 0 ; i < elongations.size() ; i++ ){
            int j = int(v_nb*(elongations[i]-min_value)/(max_value-min_value));
            values[j] ++;
        }

    }

    int VoxelGrid::getIndice(int i, int j , int k){
        return i*(zdim()+1)*(ydim()+1) + j*(zdim()+1) +k;
    }

    //void VoxelGrid::getWorldCoordinates(int i, int j, int k, BasicPoint & p){
    //    return i*(depth+1)*(height+1) + j*(depth+1) +k;
    //}

    void VoxelGrid::checkForOutliers(){

        outliers.clear();
        for(unsigned int i = 0 ; i < voxels.size() ; i ++){
            BasicPoint diff = positions[i] - getWorldCoordinate(voxels[i]);
            if(diff.sqrnorm() > 0.1){
                outliers.push_back(i);
            }
        }

        std::cout << outliers.size() << " outliers over " << voxels.size() << std::endl ;


    }

    bool VoxelGrid::ignoreOutliers(){
        /*
    checkForOutliers();
    bool found = false;
    if( outliers.size() > 0 ){
        found = true;
        
        for(unsigned int i = 0 ; i < outliers.size() ; i ++){
            setValue(voxels[outliers[i]], -1);
        }
        
        clear();
        
        for(unsigned int i = 0 ; i < _dim[0] ; i++){
            for(unsigned int j = 0 ; j < _dim[1] ; j++){
                for(unsigned int k = 0 ; k < _dim[2] ; k++){
                    GRID_TYPE label = value(i,j,k);
                    if(label >= BACKGROUND_GRID_VALUE){
                        sortedVoxels[label].push_back(voxels.size());
                        voxels.push_back(Voxel(i,j,k));
                    }
                }
            }
        }
        
        subdomain_indices.clear();
        for(std::map<int, std::vector<int> >::iterator it = sortedVoxels.begin(); it != sortedVoxels.end() ; it ++  ){
            subdomain_indices.push_back(it->first);
        }
        
        computeDefaultColors();
        
        computePositions();
        
        BasicPoint max (dx()/2., dy()/2., dz()/2.);
        BasicPoint min (-dx()/2., -dy()/2., -dz()/2.);
        
        //        initTextures( min, max );
        //        computeTranslationData();
        
        outliers.clear();
    }
    
    return found;
    */
        return false;
    }




    bool VoxelGrid::isVisiblePoint( const BasicPoint & pos ){

        float vis = dot( clippingNormal, pos - pointOnClipping );

        float cutVis[3];
        for( int i = 0 ; i < 3 ; i ++ )
            cutVis[i] = (pos[i] - cut[i])*cutDirection[i];

        if( vis < 0. || cutVis[0] < 0.|| cutVis[1] < 0.|| cutVis[2] < 0. )
            return false;
        return true;
    }

