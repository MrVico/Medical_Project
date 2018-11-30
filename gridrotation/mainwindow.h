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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "VoxelGrid.h"
#include "TetMeshCreator.h"
#include <QMainWindow>
#include <QDir>
#define DEFAULT_THRESHOLD 5
#define MAX_THRESHOLD 254
#define DMAX 10.
#define DSTEP 0.05


template <typename T>
struct GridAdjustOptionsStruct {
    T dx;
    T dy;
    T dz;
    int threshold;
    unsigned int width;
    unsigned int height;
    unsigned int depth;
    T wwidth;
    T wheight;
    T wdepth;
    bool fitToData;
    int interpolation; // 0 Nearest neighbor, 1 Trilinear
};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();
    virtual ~MainWindow();
    VoxelGrid HRGrid;
    VoxelGrid HRShearedGrid;
    TetMeshCreator tetMeshCreator;
    void openSlices( const QString & foldername, GridAdjustOptionsStruct<double> &gridAdjustInputOptions );
    void saveSlices( const QString & foldername, const GridAdjustOptionsStruct<double> & gridAdjustOptions  );
    void saveToIMA( const QDir & dir, VoxelGrid & grid, int threshold);
    bool inputParameters( GridAdjustOptionsStruct<double> &gridAdjustInputOptions, bool output=false );
    void updateThresholdTMBounds(VoxelGrid & HRGrid, std::ifstream & inputImaFile );
    void computeTM();
protected:
    QString _inputfolderName;
    QString _tmp_grid_ima_path;
    QString _tmp_grid_dim_path;

    Voxel _min;
    Voxel _max;

    int _threshold;

    void clear();
private slots:
    void onAddNew();
    void openSlices();
    void saveSlices();
    void processSlices();
};

#endif
