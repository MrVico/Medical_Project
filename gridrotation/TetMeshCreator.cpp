#include "TetMeshCreator.h"
#include <QImage>
#include <QDir>
#include <QTime>
#include <fstream>


void TetMeshCreator::updateMesh(const std::vector<BasicPoint> & vertices , const std::vector<Tetrahedron> & tetrahedra, float xshearstep ){

    _vertices = vertices;
    _tetrahedra = tetrahedra;

    initial_vertices_positions = vertices;

    _xshearstep = xshearstep;
    computeBB();
}

void TetMeshCreator::updateVertices( std::vector<BasicPoint> & points ){

    _vertices = points;
    computeBB();

}

void TetMeshCreator::computeBB(){

    BBMin = BasicPoint( FLT_MAX, FLT_MAX, FLT_MAX);
    BBMax = BasicPoint(-FLT_MAX,-FLT_MAX,-FLT_MAX);

    for( unsigned int i = 0 ; i < _vertices.size() ; i ++ )
    {
        const BasicPoint & point = _vertices[i];
        for(unsigned int j = 0 ; j < 3 ; j ++ ){
            if(BBMin[j] > point[j]) BBMin[j] = point[j];
            if(BBMax[j] < point[j]) BBMax[j] = point[j];
        }
    }
}

void TetMeshCreator::computeBarycentricCoordinates(int tet_id , const BasicPoint & point, float & ld0, float & ld1, float & ld2, float & ld3 ){

    BasicPoint p1 = _vertices[_tetrahedra[tet_id].getVertex(0)];
    BasicPoint p2 = _vertices[_tetrahedra[tet_id].getVertex(1)];
    BasicPoint p3 = _vertices[_tetrahedra[tet_id].getVertex(2)];
    BasicPoint p4 = _vertices[_tetrahedra[tet_id].getVertex(3)];

    float a = p1[0] - p4[0]; float b = p2[0] - p4[0]; float c = p3[0] - p4[0];
    float d = p1[1] - p4[1]; float e = p2[1] - p4[1]; float f = p3[1] - p4[1];
    float g = p1[2] - p4[2]; float h = p2[2] - p4[2]; float k = p3[2] - p4[2];

    float det = a*(e*k - h*f) + b*(g*f - d*k) + c*(d*h - g*e);

    float A = e*k - f*h; float D = c*h - b*k; float G = b*f - c*e;
    float B = f*g - d*k; float E = a*k - c*g; float H = c*d - a*f;
    float C = d*h - e*g; float F = g*b - a*h; float K = a*e - b*d;

    ld0 = ((point[0] - p4[0])*A + (point[1] - p4[1])*D + (point[2] - p4[2])*G )/ det;
    ld1 = ((point[0] - p4[0])*B + (point[1] - p4[1])*E + (point[2] - p4[2])*H )/ det;
    ld2 = ((point[0] - p4[0])*C + (point[1] - p4[1])*F + (point[2] - p4[2])*K )/ det;
    ld3 = 1 - ld0 - ld1 - ld2;



}

void TetMeshCreator::computePseudoRasterDeformation(VoxelGrid & HRGrid, VoxelGrid & HRShearedGrid, std::ifstream & inputImaFile, VoxelGrid & HRresult, const QString & outputFoldername, int threshold, int interpolation, QStatusBar *statusbar ){

    std::cout << "Compute deformation " << HRresult.xdim() << " - " << HRresult.ydim() << " - " << HRresult.zdim() << std::endl;
    QTime t;
    t.start();

    char current_icharacter;

    // Seek to the beginning of the files
    inputImaFile.seekg(0, std::ios::beg);

    std::cout << "Compute deformation " << HRresult.xdim() << " - " << HRresult.ydim() << " - " << HRresult.zdim() << std::endl;
    std::cout << "d " << HRresult.dx() << " - " << HRresult.dy() << " - " << HRresult.dz() << std::endl;

    std::cout << "Saving slices" << std::endl;
    //  clearImage();


    QImage zero_image  ( HRresult.xdim(),HRresult.ydim(), QImage::Format_RGB32 );

    for( int i = 0 ; i < HRresult.xdim(); i ++ ){
        for( int j = 0; j < HRresult.ydim() ; j ++ ){
            zero_image.setPixel(i,j, qRgb(0, 0, 0));
        }
    }

    std::cout << HRresult.xdim() << " and " << HRresult.ydim() << " and " << HRresult.zdim() << std::endl;

    BasicPoint x6 = _vertices[5];
    BasicPoint x1 = _vertices[0];
    //Reducing the computation, add min translation ?
    float step = (x6[0])/HRresult.zdim();

    Voxel vmax = HRresult.getGridCoordinate(x1);

    int scaled_z = HRGrid.zdim();

    double d[3]= {(double)0.65,(double)0.65,(double)0.65};
    float z_shear = _xshearstep/HRGrid.dx();
    // if(result.isInGrid(vmin) && result.isInGrid(vmax)){
    float maxshear = scaled_z*z_shear - 1;
    // if(result.isInGrid(vmin) && result.isInGrid(vmax)){
    for( int k = 0 ; k < HRresult.zdim(); k ++ ){

        std::cout << "Processing image " << k << std::endl;
        QString message = QString::fromLatin1("Processing image %1/%2").arg(k).arg(HRresult.zdim());

        statusbar->showMessage(message);

        QImage current_z_image = zero_image;

        BasicPoint pmin (step*k, 0.,0.);

        Voxel vmin = HRresult.getGridCoordinate(pmin);

        for( int j = 0; j < HRresult.ydim() ; j ++ ){
//            for( int i = 0; i < HRresult.xdim() ; i ++ ){
                for( int i = vmin.i(); i <= vmin.i()+vmax.i() ; i ++ ){

                //   current_z_image.setPixel((int)i,(int)j, qRgb(255, 0, 0));

                BasicPoint position = HRresult.getWorldCoordinate(i,j,k);
                // std::cout << position << std::endl;
                bool found = false;
                int tet_id = 0 ;
                float ld[4];

                while( !found && tet_id < _tetrahedra.size() ){
                    computeBarycentricCoordinates(tet_id, position, ld[0], ld[1], ld[2], ld[3]);
                    //std::cout << ld[0] << " "<< ld[1] << " "<<ld[2] << " "<<ld[3] << std::endl;
                    if(ld[0] >= 0 && ld[0] <=1. && ld[1] >= 0 && ld[1] <=1. && ld[2] >= 0 && ld[2] <=1. && ld[3] >= 0 && ld[3] <=1. ){
                        found = true;
                    } else
                        tet_id++;
                }

                if( found ){

                    BasicPoint p0 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(0)];
                    BasicPoint p1 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(1)];
                    BasicPoint p2 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(2)];
                    BasicPoint p3 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(3)];

                    BasicPoint pos = ld[0]*p0 + ld[1]*p1 + ld[2]*p2 + ld[3]*p3;

#if 1

                    Voxel projected_v= HRShearedGrid.getGridCoordinate(pos);

                    if( projected_v.k() < HRShearedGrid.zdim() ){

                        if( projected_v.i() < HRShearedGrid.xdim() && projected_v.i() < HRShearedGrid.ydim() ){

                            int voxel_value = 0;
                            int projected_voxel_value= getValueFromFile( projected_v.i(), projected_v.j(), projected_v.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);

                            if( interpolation == 1 ){

                                BasicPoint center_coord = HRShearedGrid.getWorldCoordinate(projected_v);
                                BasicPoint vec = pos-center_coord;

                                if(vec[0] < 0.000001 & vec[1] < 0.000001 & vec[2] < 0.000001 ){
                                    voxel_value= projected_voxel_value;
                                } else {

                                    int xmin = projected_v.i();
                                    int ymin = projected_v.j();
                                    int zmin = projected_v.k();


                                    int xmax = projected_v.i();
                                    int ymax = projected_v.j();
                                    int zmax = projected_v.k();



                                    if( vec[0] <= 0 ) xmin = projected_v.i() - 1;
                                    else xmax = projected_v.i() + 1;

                                    if( vec[1] <= 0 ) ymin = projected_v.j() - 1;
                                    else ymax = projected_v.j() + 1;

                                    if( vec[2] <= 0 ) zmin = projected_v.k() - 1;
                                    else zmax = projected_v.k() + 1;

                                    Voxel v_born_min (xmin, ymin, zmin);
                                    Voxel v_born_max (xmax, ymax, zmax);

                                    BasicPoint p_0 = HRShearedGrid.getWorldCoordinate(v_born_min);
                                    BasicPoint p_1 = HRShearedGrid.getWorldCoordinate(v_born_max);

                                    float x_d = (pos[0]-p_0[0])/(p_1[0]-p_0[0]);
                                    float y_d = (pos[1]-p_0[1])/(p_1[1]-p_0[1]);
                                    float z_d = (pos[2]-p_0[2])/(p_1[2]-p_0[2]);

                                    //                                   std::cout << x_d << " , " <<  y_d <<", " << z_d << std::endl << std::endl;

                                    int c000 = getValueFromFile( v_born_min.i(), v_born_min.j(), v_born_min.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);
                                    int c100 = getValueFromFile( v_born_max.i(), v_born_min.j(), v_born_min.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);
                                    int c110 = getValueFromFile( v_born_max.i(), v_born_max.j(), v_born_min.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);
                                    int c101 = getValueFromFile( v_born_max.i(), v_born_min.j(), v_born_max.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);
                                    int c010 = getValueFromFile( v_born_min.i(), v_born_max.j(), v_born_min.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);
                                    int c011 = getValueFromFile( v_born_min.i(), v_born_max.j(), v_born_max.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);
                                    int c001 = getValueFromFile( v_born_min.i(), v_born_min.j(), v_born_max.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);
                                    int c111 = getValueFromFile( v_born_max.i(), v_born_max.j(), v_born_max.k(), HRGrid, HRShearedGrid, inputImaFile, maxshear, z_shear);



                                    if(v_born_min.i() < 0  || v_born_min.i() == HRShearedGrid.xdim() ){
                                        c000 = projected_voxel_value;
                                        c010 = projected_voxel_value;
                                        c011 = projected_voxel_value;
                                        c001 = projected_voxel_value;
                                    }

                                    if(v_born_min.j() < 0  || v_born_min.j() == HRShearedGrid.xdim() ){
                                        c000 = projected_voxel_value;
                                        c101 = projected_voxel_value;
                                        c100 = projected_voxel_value;
                                        c001 = projected_voxel_value;
                                    }

                                    if(v_born_min.k() < 0  || v_born_min.k() == HRShearedGrid.xdim() ){
                                        c000 = projected_voxel_value;
                                        c110 = projected_voxel_value;
                                        c100 = projected_voxel_value;
                                        c010 = projected_voxel_value;
                                    }

                                    if(v_born_max.i() < 0  || v_born_max.i() == HRShearedGrid.xdim() ){
                                        c100 = projected_voxel_value;
                                        c110 = projected_voxel_value;
                                        c111 = projected_voxel_value;
                                        c101 = projected_voxel_value;
                                    }

                                    if(v_born_max.j() < 0  || v_born_max.j() == HRShearedGrid.xdim() ){
                                        c010 = projected_voxel_value;
                                        c111 = projected_voxel_value;
                                        c110 = projected_voxel_value;
                                        c011 = projected_voxel_value;
                                    }

                                    if(v_born_max.k() < 0  || v_born_max.k() == HRShearedGrid.xdim() ){
                                        c001 = projected_voxel_value;
                                        c111 = projected_voxel_value;
                                        c101 = projected_voxel_value;
                                        c011 = projected_voxel_value;
                                    }

                                    // std::cout << c000 << " , " <<  c100 <<", " << c110 << ", " <<c101 << " , " <<  c010 <<", " << c011 << ", "<< c001 << ", "<< c111 << std::endl << std::endl;

                                    float c00 = c000*(1.-x_d)+c100*x_d;
                                    float c01 = c001*(1.-x_d)+c101*x_d;
                                    float c10 = c010*(1.-x_d)+c110*x_d;
                                    float c11 = c011*(1.-x_d)+c111*x_d;

                                    float c0 = c00*(1.-y_d)+c10*y_d;
                                    float c1 = c01*(1.-y_d)+c11*y_d;

                                    voxel_value = int( c0*(1.-z_d)+c1*z_d );
                                    //   voxel_value = int( float(c000 + c100 + c110 + c101+c010+c011+c001+c111)/8. );
                                    //  std::cout << voxel_value << std::endl;

                                    //  voxel_value = c100;
                                }

                            } else {

                                voxel_value = projected_voxel_value;
                            }

                            if( voxel_value > threshold ){

                                current_z_image.setPixel((int)i,(int)j, qRgb(voxel_value, voxel_value, voxel_value));
                            }
                        }

                    }
#endif
                }

            }
        }

        //}

        int num = k;

        QString count;
        QString fname = outputFoldername;
        fname.append("/");
        if(num < 10) fname.append("000");
        else if(num < 100) fname.append("00");
        else if(num < 1000) fname.append("0");

        count.setNum(num);
        fname.append(count);

        fname.append(".tif");
        current_z_image.save(fname);
        std::cout << k << " saved in " << fname.toStdString() << std::endl;
    }
    std::cout << "Time elapsed for rasterization: " << t.elapsed() << " ms" << std::endl;
}

int TetMeshCreator::getValueFromFile( int projected_v_i, int projected_v_j, int projected_v_k, VoxelGrid & HRGrid, VoxelGrid & HRShearedGrid, std::ifstream & inputImaFile, float maxshear, float z_shear){


    if(projected_v_i < 0 || projected_v_j <0 || projected_v_k < 0 || projected_v_i  == HRShearedGrid.xdim() || projected_v_j== HRShearedGrid.ydim()|| projected_v_k == HRShearedGrid.zdim()){
        return 0;
    }

    int voxel_value = 0;

    float step_f = maxshear - z_shear*projected_v_k;
    int step = (int)step_f;

    int image_i = projected_v_i - step;
    if( image_i < HRGrid.xdim()  ){
        char current_icharacter;
        int input_index = HRGrid.getGridIndex( image_i, projected_v_j, projected_v_k );
        inputImaFile.seekg(input_index, std::ios::beg);
        inputImaFile.read( &current_icharacter, 1 );
        voxel_value = (int)current_icharacter;
    }
    return voxel_value;
}

void TetMeshCreator::computePseudoRasterDeformation(VoxelGrid & HRGrid, VoxelGrid & HRShearedGrid, std::ifstream & inputImaFile, VoxelGrid & HRresult, std::ofstream & outputImaFile, int threshold ){

    std::cout << "Compute deformation " << HRresult.xdim() << " - " << HRresult.ydim() << " - " << HRresult.zdim() << std::endl;
    QTime t;
    t.start();

    char current_icharacter;

    // Seek to the beginning of the files
    inputImaFile.seekg(0, std::ios::beg);

    std::cout << HRresult.xdim() << " and " << HRresult.ydim() << " and " << HRresult.zdim() << std::endl;

    BasicPoint x6 = _vertices[5];
    BasicPoint x1 = _vertices[0];
    float step = x6[0]/HRresult.zdim();

    Voxel vmax = HRresult.getGridCoordinate(x1);

    // if(result.isInGrid(vmin) && result.isInGrid(vmax)){

    // if(result.isInGrid(vmin) && result.isInGrid(vmax)){
    for( int k = 0 ; k < HRresult.zdim(); k ++ ){

        std::cout << "Processing image " << k << std::endl;

        BasicPoint pmin (step*k, 0.,0.);

        Voxel vmin = HRresult.getGridCoordinate(pmin);

        for( int j = 0; j < HRresult.ydim() ; j ++ ){
            for( int i = 0; i < HRresult.xdim() ; i ++ ){
                //for( int i = vmin.i(); i <= vmin.i()+vmax.i() ; i ++ ){

                //   current_z_image.setPixel((int)i,(int)j, qRgb(255, 0, 0));

                BasicPoint position = HRresult.getWorldCoordinate(i,j,k);
                // std::cout << position << std::endl;
                bool found = false;
                int tet_id = 0 ;
                float ld[4];

                while( !found && tet_id < _tetrahedra.size() ){
                    computeBarycentricCoordinates(tet_id, position, ld[0], ld[1], ld[2], ld[3]);
                    //std::cout << ld[0] << " "<< ld[1] << " "<<ld[2] << " "<<ld[3] << std::endl;
                    if(ld[0] >= 0 && ld[0] <=1. && ld[1] >= 0 && ld[1] <=1. && ld[2] >= 0 && ld[2] <=1. && ld[3] >= 0 && ld[3] <=1. ){
                        found = true;
                    } else
                        tet_id++;
                }

                if( found ){

                    BasicPoint p0 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(0)];
                    BasicPoint p1 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(1)];
                    BasicPoint p2 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(2)];
                    BasicPoint p3 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(3)];

                    BasicPoint pos = ld[0]*p0 + ld[1]*p1 + ld[2]*p2 + ld[3]*p3;

                    // std::cout << pos << std::endl;

                    Voxel projected_v= HRGrid.getGridCoordinate(pos);

                    if( projected_v.k() < HRGrid.zdim() ){

                        if( projected_v.i() < HRGrid.xdim() && projected_v.i() < HRGrid.ydim() ){
                            int input_index = HRGrid.getGridIndex( projected_v );

                            // std::cout << pos << std::endl;


                            inputImaFile.seekg(input_index, std::ios::beg);
                            inputImaFile.read( &current_icharacter, 1 );

                            if( (int)current_icharacter > threshold ){

                                int output_index = HRresult.getGridIndex(i,j,k);

                                outputImaFile.seekp(output_index, std::ios::beg);
                                outputImaFile.write( &current_icharacter, 1 );
                            }
                            // si = 1;


                        }
                    }

                }
            }
        }
    }
    //}

    std::cout << "Time elapsed for rasterization: " << t.elapsed() << " ms" << std::endl;
    // result.sort();

    //updatePositions();

}



void TetMeshCreator::computePseudoRasterDeformation(VoxelGrid & HRGrid, VoxelGrid & HRShearedGrid, QString & foldername, VoxelGrid & HRresult, const QString & outputFoldername, int threshold ){

    std::cout << "Compute deformation " << HRresult.xdim() << " - " << HRresult.ydim() << " - " << HRresult.zdim() << std::endl;
    std::cout << "d " << HRresult.dx() << " - " << HRresult.dy() << " - " << HRresult.dz() << std::endl;
    QTime t;
    t.start();

    QDir dir (foldername);

    dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
    dir.setSorting(QDir::Size | QDir::Reversed);


    QStringList filters;
    filters << "*.png" ;
    filters << "*.tif" ;
    dir.setNameFilters(filters);
    dir.setSorting(QDir::Name);
    QFileInfoList list = dir.entryInfoList();

    if( list.size() < 1 ) return ;

    std::cout << "Coucou open slices!!!!!!" << std::endl;
    //  clearImage();

    QFileInfo fileInfo = list.at(0);
    QImage zero_image  ( fileInfo.absoluteFilePath() );

    zero_image = zero_image.scaled(HRresult.xdim(),HRresult.ydim(),Qt::IgnoreAspectRatio,Qt::SmoothTransformation);

    for( int i = 0 ; i < HRresult.xdim(); i ++ ){
        for( int j = 0; j < HRresult.ydim() ; j ++ ){
            zero_image.setPixel(i,j, qRgb(0, 0, 0));
        }
    }

    std::cout << HRresult.xdim() << " and " << HRresult.ydim() << " and " << HRresult.zdim() << std::endl;

    BasicPoint x6 = _vertices[5];
    BasicPoint x1 = _vertices[0];
    float step = x6[0]/HRresult.zdim();

    std::vector<QImage> images;
    for (int k = 0; k < list.size(); ++k) {
        QFileInfo fileInfo = list.at(k);
        QImage img  ( fileInfo.absoluteFilePath() );
        images.push_back(img);
    }


    Voxel vmax = HRresult.getGridCoordinate(x1);

    // if(result.isInGrid(vmin) && result.isInGrid(vmax)){

    // if(result.isInGrid(vmin) && result.isInGrid(vmax)){
    for( int k = 0 ; k < HRresult.zdim(); k ++ ){

        std::cout << "Processing image " << k << std::endl;
        QImage current_z_image = zero_image;

        BasicPoint pmin (step*k, 0.,0.);

        Voxel vmin = HRresult.getGridCoordinate(pmin);

        for( int j = 0; j < HRresult.ydim() ; j ++ ){
            for( int i = 0; i < HRresult.xdim() ; i ++ ){
                //for( int i = vmin.i(); i <= vmin.i()+vmax.i() ; i ++ ){

                //   current_z_image.setPixel((int)i,(int)j, qRgb(255, 0, 0));

                BasicPoint position = HRresult.getWorldCoordinate(i,j,k);
                position += BasicPoint(0.000001, 0.000001, 0.000001);
                // std::cout << position << std::endl;
                bool found = false;
                int tet_id = 0 ;
                float ld[4];

                while( !found && tet_id < _tetrahedra.size() ){
                    computeBarycentricCoordinates(tet_id, position, ld[0], ld[1], ld[2], ld[3]);
                    //std::cout << ld[0] << " "<< ld[1] << " "<<ld[2] << " "<<ld[3] << std::endl;
                    if(ld[0] >= 0 && ld[0] <=1. && ld[1] >= 0 && ld[1] <=1. && ld[2] >= 0 && ld[2] <=1. && ld[3] >= 0 && ld[3] <=1. ){
                        found = true;
                    } else
                        tet_id++;
                }

                if( found ){

                    BasicPoint p0 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(0)];
                    BasicPoint p1 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(1)];
                    BasicPoint p2 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(2)];
                    BasicPoint p3 = initial_vertices_positions[_tetrahedra[tet_id].getVertex(3)];

                    BasicPoint pos = ld[0]*p0 + ld[1]*p1 + ld[2]*p2 + ld[3]*p3;

                    // std::cout << pos << std::endl;

                    Voxel projected_v= HRGrid.getGridCoordinate(pos);

                    if( projected_v.k() < HRGrid.zdim() ){

                        QFileInfo fileInfo = list.at(projected_v.k());
                        QImage projected_img  ( fileInfo.absoluteFilePath() );

                        if( projected_v.i() < HRGrid.xdim() && projected_v.i() < HRGrid.ydim() ){
                            QColor color = projected_img.pixel( projected_v.i(), projected_v.j() );

                            if( color.red() >  threshold ){
                                current_z_image.setPixel((int)i,(int)j, color.rgb());

                            }
                        }
                    }
                }

            }
        }

        //}

        int num = k;

        QString count;
        QString fname = outputFoldername;
        fname.append("/");
        if(num < 10) fname.append("000");
        else if(num < 100) fname.append("00");
        else if(num < 1000) fname.append("0");

        count.setNum(num);
        fname.append(count);

        fname.append(".tif");
        current_z_image.save(fname);
        std::cout << k << " saved in " << fname.toStdString() << std::endl;
    }
    std::cout << "Time elapsed for current rasterization: " << t.elapsed() << " ms" << std::endl;
}
