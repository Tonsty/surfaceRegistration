#include "PCAStatistics.h"

double ** PCAStatistics(vtkSmartPointer<vtkPolyData> data, int pointCount, int orientation) {
    
    // 1. Get some data
    
    double rowDataAdjust[3][pointCount];
    
    const char m0Name[] = "x";
    vtkSmartPointer<vtkDoubleArray> dataSetX = vtkSmartPointer<vtkDoubleArray>::New();
    dataSetX->SetNumberOfComponents(1);
    dataSetX->SetName( m0Name );
    
    const char m1Name[] = "y";
    vtkSmartPointer<vtkDoubleArray> dataSetY = vtkSmartPointer<vtkDoubleArray>::New();
    dataSetY->SetNumberOfComponents(1);
    dataSetY->SetName( m1Name );
    
    const char m2Name[] = "z";
    vtkSmartPointer<vtkDoubleArray> dataSetZ = vtkSmartPointer<vtkDoubleArray>::New();
    dataSetZ->SetNumberOfComponents(1);
    dataSetZ->SetName( m2Name );
    
    // 2. Subtract the mean
    
    double countX = 0, countY = 0, countZ = 0;
    
    for(vtkIdType i = 0; i < pointCount; i++)
    {
        double p[3];
        data->GetPoint(i,p);
        countX += p[0];
        countY += p[1];
        countZ += p[2];
    }
    
    double meanX = countX / pointCount;
    double meanY = countY / pointCount;
    double meanZ = countZ / pointCount;
    
    // Create the RowDataAdjust
    
    for(vtkIdType i = 0; i < pointCount; i++)
    {
        double p[3];
        data->GetPoint(i,p);
        dataSetX->InsertNextValue(p[0]);
        dataSetY->InsertNextValue(p[1]);
        dataSetZ->InsertNextValue(p[2]);
        
        rowDataAdjust[0][i] = p[0] - meanX;
        rowDataAdjust[1][i] = p[1] - meanY;
        rowDataAdjust[2][i] = p[2] - meanZ;
    }
    
    vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
    datasetTable->AddColumn(dataSetX);
    datasetTable->AddColumn(dataSetY);
    datasetTable->AddColumn(dataSetZ);
    
    // 4. Calculate the eigenvectors and eigenvalues of the covariance matrix
    
    vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
#if VTK_MAJOR_VERSION <= 5
    pcaStatistics->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable );
#else
    pcaStatistics->SetInputData( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable );
#endif
    
    pcaStatistics->SetColumnStatus("x", 1 );
    pcaStatistics->SetColumnStatus("y", 1 );
    pcaStatistics->SetColumnStatus("z", 1 );
    pcaStatistics->RequestSelectedColumns();
    pcaStatistics->SetDeriveOption(true);
    pcaStatistics->Update();
    
    vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
    pcaStatistics->GetEigenvectors(eigenvectors);
    
    vtkSmartPointer<vtkDoubleArray> evec1 = vtkSmartPointer<vtkDoubleArray>::New();
    pcaStatistics->GetEigenvector(0, evec1);
    
    vtkSmartPointer<vtkDoubleArray> evec2 = vtkSmartPointer<vtkDoubleArray>::New();
    pcaStatistics->GetEigenvector(1, evec2);
    
    vtkSmartPointer<vtkDoubleArray> evec3 = vtkSmartPointer<vtkDoubleArray>::New();
    pcaStatistics->GetEigenvector(2, evec3);
    
    // 5. Choosing components and forming a feature vector
    
    double rowFeatureVector[3][3];
    
    for (int i = 0; i < 3; i++) {
        rowFeatureVector[0][i] = (orientation == 1 || orientation == 4 || orientation == 5 || orientation == 7) ? -evec1->GetValue(i) : evec1->GetValue(i);
        rowFeatureVector[1][i] = (orientation == 2 || orientation == 4 || orientation == 6 || orientation == 7) ? -evec2->GetValue(i) : evec2->GetValue(i);
        rowFeatureVector[2][i] = (orientation == 3 || orientation == 5 || orientation == 6 || orientation == 7) ? -evec3->GetValue(i) : evec3->GetValue(i);
    }
    
    // 6. Deriving the dataset
    
    double **finalData = new double*[3];
    
    for (int i = 0; i != 3; ++i)
    {
        finalData[i] = new double[pointCount];
        for (int j = 0; j != pointCount; ++j)
        {
            double sum = 0;
            for (int k = 0; k != 3; ++k)
            {
                sum += rowFeatureVector[i][k] * rowDataAdjust[k][j];
            }
            finalData[i][j] = sum;
        }
    }
    
    return finalData;
}

