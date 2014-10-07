#include "preProcessing.h"

void freeData(double **finalData)
{
    for (int i=0; i<3; i++)
        delete [] finalData[i];
    delete [] finalData;
}

vtkSmartPointer<vtkPolyData> TransformScaleTranslate(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkPolyData> polydata, double value)
{
    vtkSmartPointer<vtkTransform> customTransform = vtkSmartPointer<vtkTransform>::New();
    customTransform->Scale(value, value, value);
    
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
#if VTK_MAJOR_VERSION <= 5
    transformFilter->SetInput(polydata);
#else
    transformFilter->SetInputData(polydata);
#endif
    transformFilter->SetTransform(customTransform);
    transformFilter->Update();
    
    return transformFilter->GetOutput();
}

vtkSmartPointer<vtkOctreePointLocator> buildOctree(vtkSmartPointer<vtkPolyData> newdata, vtkSmartPointer<vtkPoints> points,
                                                   vtkSmartPointer<vtkPolyData> olddata) {
    newdata->SetPoints(points);
    newdata->SetPolys(olddata->GetPolys());
    
    vtkSmartPointer<vtkOctreePointLocator> octree = vtkSmartPointer<vtkOctreePointLocator>::New();
    octree->SetDataSet(newdata);
    octree->BuildLocator();
    
    return octree;
}

double preProcessingStep(int pointCountTarget, int pointCountTemplate,
                         vtkSmartPointer<vtkPolyData> dataTarget, vtkSmartPointer<vtkPolyData> dataTemplate,
                         vtkSmartPointer<vtkPolyData> polydataTarget, vtkSmartPointer<vtkPolyData> polydataTemplate, vtkSmartPointer<vtkPolyData> transformedTemplate,
                         int maxDistance) {
    
    double maxCoordinate = 0;
    
    double **finalDataTarget;
    finalDataTarget = PCAStatistics(dataTarget, pointCountTarget, 0);
    
    vtkSmartPointer<vtkPoints> pointsTarget = vtkSmartPointer<vtkPoints>::New();
    
    for(int i = 0; i < pointCountTarget; ++i)
    {
        pointsTarget->InsertNextPoint ( finalDataTarget[0][i], finalDataTarget[1][i], finalDataTarget[2][i] );
        for (int j=0; j<3; j++) {
            if (abs(finalDataTarget[j][i]) > maxCoordinate) {
                maxCoordinate = abs(finalDataTarget[j][i]);
            }
        }
    }
    
    // Building the octree
    
    vtkSmartPointer<vtkOctreePointLocator> octreeTarget = buildOctree(polydataTarget, pointsTarget, dataTarget);
    
    // Used to rotate the template after the PCA if needed
    
    double **finalDataTemplate;
    double distance;
    double overallDistance = 0;
    double minDistance = 1000000;
    int minOrientation = 0;
    
    vtkSmartPointer<vtkPoints> pointsTemplate = vtkSmartPointer<vtkPoints>::New();
    
    for (int orientation = 0; orientation < 8; orientation++) {
        
        finalDataTemplate = PCAStatistics(dataTemplate, pointCountTemplate, orientation);
        
        for(int i = 0; i < pointCountTemplate; ++i)
        {
            pointsTemplate->InsertNextPoint( finalDataTemplate[0][i], finalDataTemplate[1][i], finalDataTemplate[2][i] );
        }
        
        for(int i = 0; i < pointCountTemplate; ++i)
        {
            octreeTarget->FindClosestPoint( finalDataTemplate[0][i], finalDataTemplate[1][i], finalDataTemplate[2][i], distance );
            overallDistance += distance;
        }
        
        if (minDistance > overallDistance) {
            minDistance = overallDistance;
            minOrientation = orientation;
        }
        
        overallDistance = 0;
        pointsTemplate = vtkSmartPointer<vtkPoints>::New();
    }
    
    finalDataTemplate = PCAStatistics(dataTemplate, pointCountTemplate, minOrientation);
    
    // The initialization of the transformed template is the same as the actual template
    
    vtkSmartPointer<vtkPoints> transformedPointsTemplate = vtkSmartPointer<vtkPoints>::New();
    
    for(int i = 0; i < pointCountTemplate; ++i)
    {
        pointsTemplate->InsertPoint( i, finalDataTemplate[0][i], finalDataTemplate[1][i], finalDataTemplate[2][i] );
        transformedPointsTemplate->InsertPoint( i, finalDataTemplate[0][i], finalDataTemplate[1][i], finalDataTemplate[2][i] );
        for (int j=0; j<3; j++) {
            if (abs(finalDataTemplate[j][i]) > maxCoordinate) {
                maxCoordinate = abs(finalDataTemplate[j][i]);
            }
        }
    }
    
    freeData(finalDataTarget);
    freeData(finalDataTemplate);
    
    // Initialization of all polydata
    
    polydataTemplate->SetPoints(pointsTemplate);
    polydataTemplate->SetPolys(dataTemplate->GetPolys());
    
    transformedTemplate->SetPoints(transformedPointsTemplate);
    transformedTemplate->SetPolys(dataTemplate->GetPolys());
    
    return maxCoordinate;
}