/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "mainwindow.h"
#include "window.h"
#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <fstream>
#include <QFileDialog>
#include <QStatusBar>
#include <iostream>
#include <QVBoxLayout>
#include <QLabel>
#include <QGroupBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QComboBox>

#define _USE_MATH_DEFINES
#include <cmath>

MainWindow::MainWindow()
{
    QMenuBar *menuBar = new QMenuBar;
    QMenu *menuWindow = menuBar->addMenu(tr("&Window"));
    QAction * fileOpen3DImageAction = new QAction (QPixmap ("./Icons/fileopenimage.png"), "Open 3D Image", this);
    fileOpen3DImageAction->setShortcut (tr ("Ctrl+I"));
    connect (fileOpen3DImageAction, SIGNAL (triggered ()) , this, SLOT (openSlices()));
    menuWindow->addAction(fileOpen3DImageAction);

    QAction * fileSaveGridAction = new QAction (QPixmap ("./Icons/filesaveimage.png"), "Save Grid", this);
    fileSaveGridAction->setShortcut (tr ("Ctrl+Alt+S"));
    connect (fileSaveGridAction, SIGNAL (triggered ()) , this, SLOT (saveSlices ()));
    menuWindow->addAction(fileSaveGridAction);
    setMenuBar(menuBar);

    _threshold = 0;
    onAddNew();

}

MainWindow::~MainWindow(){
    if(!_tmp_grid_ima_path.isEmpty()){
        remove(_tmp_grid_dim_path.toStdString().c_str());
        remove(_tmp_grid_ima_path.toStdString().c_str());
    }
}

void MainWindow::onAddNew()
{
    if (!centralWidget())
        setCentralWidget(new Window(this));
    else
        QMessageBox::information(0, tr("Cannot add new window"), tr("Already occupied. Undock first."));
}

void MainWindow::saveSlices(){

    QString folderName = QFileDialog::getExistingDirectory(this, tr("Save all images in selected directory"),
                                                           QDir::currentPath(),
                                                           QFileDialog::ShowDirsOnly
                                                           | QFileDialog::DontResolveSymlinks);

    // In case of Cancel
    if ( folderName.isEmpty() ) {
        return;
    }

    statusBar()->showMessage("Saving slices...");
    GridAdjustOptionsStruct<double> gridAdjustOptions;
    inputParameters(gridAdjustOptions, true);

    saveSlices(folderName, gridAdjustOptions);

    statusBar()->showMessage("3D image saved");
}

void MainWindow::saveSlices( const QString & foldername, const GridAdjustOptionsStruct<double> & gridAdjustOptions ){

    std::ifstream imaFile (_tmp_grid_ima_path.toStdString().c_str(), std::ios::in|std::ios::binary|std::ios::ate);
    if (!imaFile.is_open()){
        std::cout << "Error opening the file" << std::endl;
        return;
    }

    if( gridAdjustOptions.fitToData && gridAdjustOptions.threshold != _threshold ){
        _threshold = gridAdjustOptions.threshold;
        updateThresholdTMBounds(HRGrid, imaFile );
        computeTM();
    }

    BasicPoint point = tetMeshCreator.getBBMax() - tetMeshCreator.getBBMin() ;

    int n[3];

    double d[3] = {gridAdjustOptions.dx, gridAdjustOptions.dy, gridAdjustOptions.dz};

    for(int i = 0 ; i < 3 ; i++){
        n[i] = int(fabs(point[i])/d[i])+1;
    }

    VoxelGrid HRDefGrid = VoxelGrid(n[0], n[1], n[2], d[0], d[1], d[2], false, tetMeshCreator.getBBMin());


#if 1

    tetMeshCreator.computePseudoRasterDeformation(HRGrid, HRShearedGrid, imaFile, HRDefGrid, foldername, gridAdjustOptions.threshold, gridAdjustOptions.interpolation, statusBar() );

    imaFile.close();
#else
    std::ofstream outputImaFile ("/home/noura/Documents/Data/slices/sheared/out_test.ima", std::ios::out | std::ios::binary);
    if (!outputImaFile.is_open())
        return ;

    std::ofstream dimFile ("/home/noura/Documents/Data/slices/sheared/out_test.dim");
    if (!dimFile.is_open())
        return ;

    dimFile << HRDefGrid.dim(0); dimFile << " " << HRDefGrid.dim(1); dimFile<< " "  << HRDefGrid.dim(2) << std::endl;

    dimFile << "-type U8" << std::endl;

    dimFile << "-dx "; dimFile << HRDefGrid.d(0) << std::endl;
    dimFile << "-dy "; dimFile << HRDefGrid.d(1) << std::endl;
    dimFile << "-dz "; dimFile << HRDefGrid.d(2) << std::endl;

    dimFile.close();

    char buffer = 0;


    for(unsigned int i = 0 ; i < HRDefGrid.size() ; i ++){
        outputImaFile.write (&buffer, 1);
    }


    tetMeshCreator.computePseudoRasterDeformation(HRGrid, imaFile, HRDefGrid, outputImaFile);
    outputImaFile.close();
#endif
    imaFile.close();

}

void MainWindow::computeTM(){

    std::vector<BasicPoint> TM_vertices;
    std::vector<Tetrahedron> tetrahedra;

    BasicPoint pmin (_min.i()*HRGrid.dx(), _min.j()*HRGrid.dy(), _min.k()*HRGrid.dz());
    BasicPoint pmax ((_max.i()+1)*HRGrid.dx(), (_max.j()+1)*HRGrid.dy(), (_max.k()+1)*HRGrid.dz());

    TM_vertices.push_back(BasicPoint(pmin[0], pmin[1], pmin[2]));
    TM_vertices.push_back(BasicPoint(pmin[0], pmin[1], pmax[2]));
    TM_vertices.push_back(BasicPoint(pmin[0], pmax[1], pmin[2]));
    TM_vertices.push_back(BasicPoint(pmin[0], pmax[1], pmax[2]));
    TM_vertices.push_back(BasicPoint(pmax[0], pmin[1], pmin[2]));
    TM_vertices.push_back(BasicPoint(pmax[0], pmin[1], pmax[2]));
    TM_vertices.push_back(BasicPoint(pmax[0], pmax[1], pmin[2]));
    TM_vertices.push_back(BasicPoint(pmax[0], pmax[1], pmax[2]));

    tetrahedra.push_back(Tetrahedron(4,6,2,7));
    tetrahedra.push_back(Tetrahedron(1,0,3,5));
    tetrahedra.push_back(Tetrahedron(3,4,7,5));
    tetrahedra.push_back(Tetrahedron(0,4,3,5));
    tetrahedra.push_back(Tetrahedron(0,4,2,3));
    tetrahedra.push_back(Tetrahedron(2,4,7,3));



    float z_shear_step = 5./(sqrt(2.));
    //TODO make it work with threshold: different shear : depend on z of min max

    BasicPoint translation((HRGrid.dim(2) - _min.k())*z_shear_step,0.,0.);

    TM_vertices[0] += translation;
    TM_vertices[2] += translation;
    TM_vertices[4] += translation;
    TM_vertices[6] += translation;

    translation = BasicPoint ((HRGrid.dim(2) - (_max.k()+1))*z_shear_step,0.,0.);

    TM_vertices[1] += translation;
    TM_vertices[3] += translation;
    TM_vertices[5] += translation;
    TM_vertices[7] += translation;

    tetMeshCreator = TetMeshCreator( );
    tetMeshCreator.updateMesh(TM_vertices, tetrahedra, z_shear_step);

    BasicPoint min ( FLT_MAX, FLT_MAX, FLT_MAX );
    translation = -1.*TM_vertices[1];

    for(unsigned int i=0 ; i <TM_vertices.size() ; i++){
        BasicPoint p = TM_vertices[i] + translation;
        TM_vertices[i] = BasicPoint ( p[0]*cos(-M_PI/4.) + p[2]*sin(-M_PI/4.), p[1], -1.*p[0]*sin(-M_PI/4.) + p[2]*cos(-M_PI/4.) );
    }
#if 0
    QString filename = QString("./box.mesh");
    std::ofstream out (filename.toUtf8());
    if (!out)
        exit (EXIT_FAILURE);

    out << "MeshVersionFormatted " << 1 << std::endl;
    out << "Dimension " << 3 << std::endl;
    out << "Vertices\n";
    out << 8 << std::endl;

    for( unsigned int v = 0 ; v < TM_vertices.size() ; ++v )
    {
        out << (TM_vertices[v]) ;
        out <<" " << 2 <<std::endl;
    }

    out << "Triangles\n";
    out << 0 << std::endl;

    out << "Tetrahedra\n";
    out << 6 << std::endl;
    out <<5 << " " <<7<< " " << 3<< " " << 8<< " " << 1<< std::endl;
    out <<2<< " " << 1<< " " << 4<< " " << 6<< " " << 1<< std::endl;
    out <<4<< " " << 5<< " " << 8<< " " << 6<< " " << 1<< std::endl;
    out <<1<< " " << 5<< " " << 4<< " " << 6<< " " << 1<< std::endl;
    out <<1<< " " << 5<< " " << 3<< " " << 4<< " " << 1<< std::endl;
    out <<3<< " " << 5<< " " << 8<< " " << 4<< " " << 1<< std::endl;

    out << std::endl;
    out.close ();
#endif

    tetMeshCreator.updateVertices(TM_vertices);

    std::cout << tetMeshCreator.getBBMin() << std::endl;
    std::cout << tetMeshCreator.getBBMax() << std::endl;

}

void MainWindow::processSlices(){

    _inputfolderName = QFileDialog::getExistingDirectory(this, tr("Opens all images of selected directory"),
                                                         QDir::currentPath(),
                                                         QFileDialog::ShowDirsOnly
                                                         | QFileDialog::DontResolveSymlinks);

    // In case of Cancel
    if ( _inputfolderName.isEmpty() ) {
        return;
    }

    GridAdjustOptionsStruct<double> gridAdjustInputOptions;
    inputParameters(gridAdjustInputOptions);
    _threshold = gridAdjustInputOptions.threshold;


    QString folderName = QFileDialog::getExistingDirectory(this, tr("Save all images in selected directory"),
                                                           QDir::currentPath(),
                                                           QFileDialog::ShowDirsOnly
                                                           | QFileDialog::DontResolveSymlinks);

    // In case of Cancel
    if ( folderName.isEmpty() ) {
        return;
    }

    GridAdjustOptionsStruct<double> gridAdjustOptions;
    inputParameters(gridAdjustOptions, true);

    statusBar()->showMessage("Opening slices...");
    openSlices(_inputfolderName, gridAdjustInputOptions);
    statusBar()->showMessage("3D image opened");

    statusBar()->showMessage("Saving slices...");
    saveSlices(folderName, gridAdjustOptions);
    statusBar()->showMessage("3D image saved");

}

void MainWindow::openSlices(){

    _inputfolderName = QFileDialog::getExistingDirectory(this, tr("Opens all images of selected directory"),
                                                         QDir::currentPath(),
                                                         QFileDialog::ShowDirsOnly
                                                         | QFileDialog::DontResolveSymlinks);

    // In case of Cancel
    if ( _inputfolderName.isEmpty() ) {
        return;
    }

    statusBar()->showMessage("Opening slices...");

    GridAdjustOptionsStruct<double> gridAdjustInputOptions;
    inputParameters(gridAdjustInputOptions);
    _threshold = gridAdjustInputOptions.threshold;

    statusBar()->showMessage("Opening slices...");

    openSlices(_inputfolderName, gridAdjustInputOptions);

    statusBar()->showMessage("3D image opened");

}


void MainWindow::updateThresholdTMBounds(VoxelGrid & HRGrid, std::ifstream & inputImaFile ){

    std::cout << "updateThresholdTMBounds " << _threshold << std::endl;

    char current_icharacter;

    // Seek to the beginning of the files
    inputImaFile.seekg(0, std::ios::beg);

    _max = Voxel( 0, 0, 0 );
    _min = Voxel(  HRGrid.dim(0),  HRGrid.dim(1),  HRGrid.dim(2) );

    for( int k = 0 ; k < (int)HRGrid.zdim(); k ++ ){
        for( int j = 0; j < (int)HRGrid.ydim() ; j ++ ){
            for( int i = 0; i < (int)HRGrid.xdim() ; i ++ ){

                int input_index = HRGrid.getGridIndex( i,j,k );

                inputImaFile.seekg(input_index, std::ios::beg);
                inputImaFile.read( &current_icharacter, 1 );

                if( (int)current_icharacter > _threshold ){
                    _min = Voxel::min( Voxel(i,j,k), _min );
                    _max = Voxel::max( Voxel(i,j,k), _max );

                }
            }
        }
    }

}


bool MainWindow::inputParameters( GridAdjustOptionsStruct<double> &gridAdjustInputOptions, bool output ){

    float dmax = DMAX;
    float dstep = DSTEP;

    std::cout << "input parameters" << std::endl;
    QDialog * dialog = new QDialog(this);
    dialog->setWindowTitle(tr("Slices parameters"));
    dialog->setWindowFlags(Qt::Window | Qt::CustomizeWindowHint | Qt::WindowTitleHint| Qt::WindowSystemMenuHint );
    QVBoxLayout * dilogLayout = new QVBoxLayout(dialog);

    QGroupBox* gridAdjustParameterGroupBox = new QGroupBox("Images parameters");
    dilogLayout->addWidget(gridAdjustParameterGroupBox);

    QGridLayout * gridAdjustParamGroupBoxGLayout = new QGridLayout( gridAdjustParameterGroupBox );
    int current_line = 0;

    QLabel * dxLabel = new QLabel("dx");
    QDoubleSpinBox * dxSpinBox = new QDoubleSpinBox;
    dxSpinBox->setRange(dstep, dmax);
    dxSpinBox->setSingleStep(dstep);
    dxSpinBox->setValue(0.65);
    dxSpinBox->setDecimals(2);

    gridAdjustParamGroupBoxGLayout->addWidget(dxLabel, current_line,0);
    gridAdjustParamGroupBoxGLayout->addWidget(dxSpinBox, current_line++,1);

    QLabel * dyLabel = new QLabel("dy");
    QDoubleSpinBox *dySpinBox = new QDoubleSpinBox;
    dySpinBox->setRange(dstep, dmax);
    dySpinBox->setSingleStep(dstep);
    dySpinBox->setValue(0.65);
    dySpinBox->setDecimals(2);

    gridAdjustParamGroupBoxGLayout->addWidget(dyLabel, current_line,0);
    gridAdjustParamGroupBoxGLayout->addWidget(dySpinBox, current_line++,1);

    QLabel * dzLabel = new QLabel("dz");
    QDoubleSpinBox *dzSpinBox = new QDoubleSpinBox;
    dzSpinBox->setRange(dstep, dmax);
    dzSpinBox->setSingleStep(dstep);
    dzSpinBox->setValue(5./sqrt(2.));
    dzSpinBox->setDecimals(2);

    gridAdjustParamGroupBoxGLayout->addWidget(dzLabel, current_line,0);
    gridAdjustParamGroupBoxGLayout->addWidget(dzSpinBox, current_line++,1);

    QLabel * interpolationLabel = new QLabel("Interpolation");
    QComboBox * modeComboBox = new QComboBox ();
    modeComboBox->addItem ("Nearest Neighbor");
    modeComboBox->addItem ("Trilinear interpolation");

    if(output){
        gridAdjustParamGroupBoxGLayout->addWidget(interpolationLabel, current_line,0);
        gridAdjustParamGroupBoxGLayout->addWidget(modeComboBox, current_line++,1);
    }
    QSpinBox *thresholdSpinBox = new QSpinBox;
    thresholdSpinBox->setRange(0, 254);
    thresholdSpinBox->setSingleStep(1);
    thresholdSpinBox->setValue(DEFAULT_THRESHOLD);

    thresholdSpinBox->setEnabled(true);

    gridAdjustParamGroupBoxGLayout->addWidget(new QLabel("Threshold"), current_line,0);
    gridAdjustParamGroupBoxGLayout->addWidget(thresholdSpinBox, current_line++,1);

    QCheckBox *adjustToDataCheckBox = new QCheckBox("Fit result to data");
    adjustToDataCheckBox->setChecked(true);

    dilogLayout->addWidget(adjustToDataCheckBox);

    QDialogButtonBox * buttonBox = new QDialogButtonBox();
    buttonBox->setOrientation(Qt::Horizontal);
    buttonBox->setStandardButtons(QDialogButtonBox::Ok);

    QObject::connect(buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));

    dilogLayout->addWidget(buttonBox);
    int i = dialog->exec();

    gridAdjustInputOptions.dx = dxSpinBox->value();
    gridAdjustInputOptions.dy = dySpinBox->value();
    gridAdjustInputOptions.dz = dzSpinBox->value();

    gridAdjustInputOptions.threshold = thresholdSpinBox->value();
    gridAdjustInputOptions.fitToData = adjustToDataCheckBox->isChecked();

    if(output)
        gridAdjustInputOptions.interpolation = modeComboBox->currentIndex();
    if(i == QDialog::Rejected)
        return false;

    return true;
}

void  MainWindow::clear(){

}

void MainWindow::openSlices( const QString & foldername, GridAdjustOptionsStruct<double> &gridAdjustInputOptions ){
    std::cout << "openSlices was called..." << std::endl;
    QDir dir (foldername);

    dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
    dir.setSorting(QDir::Size | QDir::Reversed);


    QStringList filters;
    filters << "*.png" ;
    filters << "*.tif" ;
    dir.setNameFilters(filters);
    dir.setSorting(QDir::Name);

#if 0

    QFileInfoList old_list = dir.entryInfoList();


    for (int k = 0; k < old_list.size(); ++k) {
        QFileInfo fileInfo = old_list.at(k);
        QImage img  ( fileInfo.absoluteFilePath() );
        QString fname = fileInfo.absoluteFilePath().left(fileInfo.absoluteFilePath().lastIndexOf("_"));
        std::cout << fname.toStdString() << std::endl;

        QString rname = fname.section('_', -1); // str == "myapp"

        fname = fname.left(fname.lastIndexOf("_"));
        int num = rname.toInt();

        QString count;

        if(num < 10) fname.append("000");
        else if(num < 100) fname.append("00");
        else if(num < 1000) fname.append("0");

        count.setNum(num);
        fname.append(count);

        fname.append(".tif");
        img.save(fname);
    }
#endif
    QFileInfoList list = dir.entryInfoList();

    if( list.size() < 1 ) return ;

    QFileInfo fileInfo = list.at(0);
    QImage img  ( fileInfo.absoluteFilePath() );
    std::cout << "Opening images " << fileInfo.absoluteFilePath().toStdString() << std::endl;

    double d[3]= {gridAdjustInputOptions.dx,gridAdjustInputOptions.dy,gridAdjustInputOptions.dz};
    int n[3] = { img.width(),
                 img.height(),
                 list.size()};

    std::cout << "(nx,dx) = ( " << n[0] << " ; " << d[0] << " ) "<< std::endl;
    std::cout << "(ny,dy) = ( " << n[1] << " ; " << d[1] << " ) "<< std::endl;
    std::cout << "(nz,dz) = ( " << n[2] << " ; " << d[2] << " ) "<< std::endl;

    int z_shear = HRGrid.dim(2)*5./(sqrt(2.));

    HRGrid = VoxelGrid( img.width(), img.height(), list.size(), gridAdjustInputOptions.dx, gridAdjustInputOptions.dy, gridAdjustInputOptions.dz, false );
    HRShearedGrid = VoxelGrid( img.width()+z_shear, img.height(), list.size(), gridAdjustInputOptions.dx, gridAdjustInputOptions.dy, gridAdjustInputOptions.dz, false );

    saveToIMA(dir, HRGrid, gridAdjustInputOptions.threshold );

    if(!gridAdjustInputOptions.fitToData){
        _min = Voxel( 0, 0, 0 );
        _max = Voxel( n[0]-1, n[1]-1, n[2]-1 );
    }

    computeTM();

}



void MainWindow::saveToIMA( const QDir & dir, VoxelGrid & grid, int threshold){

    if( grid.size() == 0 ) {
        std::cout << "GridViewer::saveToIMA(QString fileName)::Empty grid" << std::endl;
        return;
    }

    QFileInfoList list = dir.entryInfoList();

    if( list.size() < 1 ) return ;

    std::cout << "Saving images to raw file" << std::endl;

    QFileInfo fileInfo = list.at(0);

    QString fileName = dir.absolutePath();
    fileName.append("/tmp_grid.dim");

    if(!fileName.endsWith(".dim"))
        fileName.append(".dim");
    _tmp_grid_dim_path = fileName;

    QString imaName = QString(fileName);
    imaName.replace(".dim", ".ima" );
    _tmp_grid_ima_path = imaName;

    std::ofstream imaFile (imaName.toUtf8(), std::ios::out | std::ios::binary);
    if (!imaFile.is_open())
        return ;

    std::ofstream dimFile (fileName.toUtf8());
    if (!dimFile.is_open())
        return ;


    dimFile << grid.dim(0); dimFile << " " << grid.dim(1); dimFile<< " "  << grid.dim(2) << std::endl;

    dimFile << "-type U8" << std::endl;

    dimFile << "-dx "; dimFile << grid.d(0) << std::endl;
    dimFile << "-dy "; dimFile << grid.d(1) << std::endl;
    dimFile << "-dz "; dimFile << grid.d(2) << std::endl;

    dimFile.close();

    char buffer = 0;

    for(unsigned int i = 0 ; i < grid.size() ; i ++){
        imaFile.write (&buffer, 1);
    }

    _max = Voxel( 0, 0, 0 );
    _min = Voxel(  grid.dim(0),  grid.dim(1),  grid.dim(2) );


    for (int k = 0; k < list.size(); ++k) {

        QString message = QString::fromLatin1("Processing image %1/%2").arg(k).arg(list.size());

        std::cout << message.toStdString() << std::endl;
        statusBar()->showMessage(message);
        fileInfo = list.at(k);
        QImage current_image = QImage ( fileInfo.absoluteFilePath() );
        for( int i = 0 ; i < current_image.width() ; i++ ){
            for( int j = 0 ; j < current_image.height() ; j++ ){
                QColor color = current_image.pixel( i, j );
                int v = color.red();
                if( v > 0 ){
                    int output_index = grid.getGridIndex(i,j,k);
                    buffer = static_cast<char>(std::max( v, 0 ) );
                    imaFile.seekp(output_index, std::ios::beg);
                    imaFile.write( &buffer, 1 );
                    if( v > threshold ){
                        _min = Voxel::min( Voxel(i,j,k), _min );
                        _max = Voxel::max( Voxel(i,j,k), _max );
                    }
                }
            }
        }
    }

    imaFile.close();
}
