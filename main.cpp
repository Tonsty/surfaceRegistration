#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>

#include "common/inputUtilities.h"
#include "common/visualizationUtilities.h"
#include "preProcessing.h"
#include "algorithm.h"

void usage()
{
    cout << "Usage:" << endl << " <-k | -p> <-a filepath> <-e filepath> <-d distanceValue> <-x maxStiffness> <-i minStiffness> <-s stiffnessStep>" << endl;
    
    cout << "  -k   Use a .vtk file (default)\n"
    << "  -p   Use a .vtp file\n"
    << "  -a   Target image file path\n"
    << "  -e   Template image file path\n"
    << "  -d   Maximum distance between the corresponding points\n"
    << "  -x   Maximum stiffness value\n"
    << "  -i   Minimum stiffness value\n"
    << "  -s   Stiffness step value\n"
    << "\n";
    
    exit(1);
}

int main( int argc, char ** argv )
{
    enum { VTK, VTP } mode = VTK;
    
    std::string targetName;
    std::string templateName;
    double maxDistance = 0;    
    double stiffness = 0;
    double minStiffness = 0;
    double stiffnessStep = 0;
    
    int c;
    extern char *optarg;
    
    if (argc < 2) {
        usage();
        return 1;
    }
    
    while ((c = getopt(argc, argv, "kpa:e:d:x:i:s:")) != -1)
    {
        switch (c)
        {
            case 'k':
            {
                mode = VTK;
                break;
            }
                
            case 'p':
            {
                mode = VTP;
                break;
            }
            case 'a':
            {
                targetName = optarg;
                break;
            }
            case 'e':
            {
                templateName = optarg;
                break;
            }
            case 'd':
            {
                maxDistance = atof(optarg);
                break;
            }
                
            case 'x':
            {
                stiffness = atof(optarg);
                break;
            }
                
            case 'i':
            {
                minStiffness = atof(optarg);
                break;
            }
                
            case 's':
            {
                stiffnessStep = atof(optarg);
                break;
            }
            default:
            {
                usage();
                return 1;
            }
        }
    }
    
    vtkSmartPointer<vtkPolyData> dataTarget = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> dataTemplate = vtkSmartPointer<vtkPolyData>::New();
    
    switch (mode)
    {
        case VTK:
        {
            dataTarget = vtkGetInput(targetName);
            dataTemplate = vtkGetInput(templateName);
            break;
        }
        case VTP:
        {
            dataTarget = vtpGetInput(targetName);
            dataTemplate = vtpGetInput(templateName);
            break;
        }
    }
    
    int pointCountTarget = dataTarget->GetNumberOfPoints();
    int pointCountTemplate = dataTemplate->GetNumberOfPoints();
    
    vtkSmartPointer<vtkPolyData> polydataTarget = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> polydataTemplate = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> transformedTemplate = vtkSmartPointer<vtkPolyData>::New();
    
    vtkSmartPointer<vtkPoints> pointsTarget = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> pointsTemplate = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> transformedPointsTemplate = vtkSmartPointer<vtkPoints>::New();
    
    double scaleFactor = applyPCA(pointCountTarget, pointCountTemplate, dataTarget, dataTemplate,
                      polydataTarget, polydataTemplate, transformedTemplate, maxDistance);
    
    // Scale to [-1,1] Cube
    
    polydataTarget = TransformScaleTranslate(polydataTarget, scaleFactor);
    transformedTemplate = TransformScaleTranslate(transformedTemplate, scaleFactor);
    polydataTemplate = TransformScaleTranslate(polydataTemplate, scaleFactor);
    
    pointsTarget = polydataTarget->GetPoints();
    pointsTemplate = polydataTemplate->GetPoints();
    transformedPointsTemplate = transformedTemplate->GetPoints();
    
    // Building the octree for the second time with the scaled mesh
    
    vtkSmartPointer<vtkPolyData> polydataResizedTargetOctree = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkOctreePointLocator> octreeResizedTarget = buildOctree(polydataResizedTargetOctree, pointsTarget, polydataTarget);
    
    vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    
    visualizeWindow(polydataTarget, polydataTemplate, transformedTemplate, window, renderer);
    
    // Main Algorithm Loop
    
    // Extract edges
    
    vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
    switch (mode)
    {
        case VTK:
        {
            extractEdges->SetInputConnection(vtkGetSource(templateName)->GetOutputPort());
            break;
        }
        case VTP:
        {
            extractEdges->SetInputConnection(vtpGetSource(templateName)->GetOutputPort());
            break;
        }
    }
    extractEdges->Update();
    vtkCellArray* linesTemplate = extractEdges->GetOutput()->GetLines();
    
    mainAlgorithmLoop(stiffness, minStiffness, stiffnessStep, maxDistance,
                      pointsTarget, transformedPointsTemplate,
                      pointCountTemplate, octreeResizedTarget, linesTemplate, window);
    
    window->Finalize();
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->SetSize(1200, 800);
    renWin->AddRenderer(renderer);
    
    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renWin);
    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
    
    interactor->Initialize();
    interactor->Start();
    
    return 0;
}

