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

enum InputModes { VTK, VTP };

int SRApplication( InputModes mode, std::string targetName, std::string templateName,
                  double maxDistance, double maxStiffness, double minStiffness, double stiffnessStep) {
    
    // Convert to BGL undirected graph
    
    BGLUndirectedGraph dataTarget;
    BGLUndirectedGraph dataTemplate;

    clock_t start;
    double duration;
    
    start = clock();
    
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
    
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    
    std::cout<< "Duration of conversion: " << duration << " seconds" << std::endl;
    
    std::cout<< "Printing 100 vertices and all edges to see if correct: " << std::endl;
    
    typedef boost::graph_traits<BGLUndirectedGraph>::vertex_descriptor vertex_t;
    typedef boost::graph_traits<BGLUndirectedGraph>::edge_iterator edge_iter;
    
    for (vertex_t v = 0; v < 100; v++)
        cout << "point " << v << ": " << dataTarget[v].x << " " << dataTarget[v].y << " " << dataTarget[v].z << endl;
    
    edge_iter ei, ei_end;

    for (tie(ei, ei_end) = edges(dataTarget); ei != ei_end; ++ei)
        cout << "(" << source(*ei, dataTarget) << "," << target(*ei, dataTarget) << ") " << endl;
    
    // Do the PCA to both target and template?
    
    // Scale the meshes?
    
    // Build Johnson's All-Pairs Shortest Paths
    
    
    
    // Use k-medoids
    
    // Run the algorithm
    
//    vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    
//    visualizeWindow(polydataTarget, polydataTemplate, transformedTemplate, window, renderer);
//    
//    window->Finalize();
//    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
//    renWin->SetSize(1200, 800);
//    renWin->AddRenderer(renderer);
//    
//    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    interactor->SetRenderWindow(renWin);
//    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
//    
//    interactor->Initialize();
//    interactor->Start();
    
    return 0;
}

int main( int argc, char ** argv )
{
    InputModes mode = VTK;
    
    std::string targetName;
    std::string templateName;
    double maxDistance = 0;
    double maxStiffness = 0;
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
                maxStiffness = atof(optarg);
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
    
    int result = SRApplication( mode, targetName, templateName, maxDistance, maxStiffness, minStiffness, stiffnessStep );
    
    return result;
}

