#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTrackball.h>
#include <vtkVersion.h>

#include <vtkDataSetReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkExtractEdges.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkOctreePointLocator.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>

#include <vtkDoubleArray.h>
#include <vtkPCAStatistics.h>
#include <vtkTable.h>
#include <vtkIdList.h>

#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;
using namespace std;

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

void freeData(double **finalData)
{
    for (int i=0; i<3; i++)
        delete [] finalData[i];
    delete [] finalData;
}

int main()
{
    
//  The source files are defined. The data is handled by vtkDataSetReader class and its output is saved in a vtkPolyData object.
//  The points and their count are extracted in vtkPoints objects. This is done for both target and template meshes.
//  Principal component analysis is done for the target mesh and an octree is build with as an object of the vtkOctreePointLocator class.
//  PCA is also done for the template mesh, with a consideration of the orientation. The final data is stored in finalDataTemplate array.
    
    int errorInt = 0;
    
    const char* targetNameInput = "/Users/sisip/Documents/Aachen University/Master Thesis/code/data/mobius.vtk";
    const char* templateNameInput = "/Users/sisip/Documents/Aachen University/Master Thesis/code/data/mobius-transformed.vtk";
    
    string targetName;
    string templateName;
    cout << "Please give the Target file path (or just press Enter to go on with the fixed example):" << endl;
    getline(cin, targetName);
    cout << "Please give the Template file path (or just press Enter to go on with the fixed example):" << endl;
    getline(cin, templateName);
    
    if (targetName.size()) {
        targetNameInput = targetName.c_str();
    }
    if (templateName.size()) {
        templateNameInput = templateName.c_str();
    }
    
    vtkSmartPointer<vtkDataSetReader> sourceTarget = vtkSmartPointer<vtkDataSetReader>::New();
    vtkSmartPointer<vtkDataSetReader> sourceTemplate = vtkSmartPointer<vtkDataSetReader>::New();
//    vtkSmartPointer<vtkXMLPolyDataReader> sourceTarget = vtkSmartPointer<vtkXMLPolyDataReader>::New();
//    vtkSmartPointer<vtkXMLPolyDataReader> sourceTemplate = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    
    sourceTarget->SetFileName(targetNameInput);
    sourceTarget->Update();
    
    sourceTemplate->SetFileName(templateNameInput);
    sourceTemplate->Update();
    
    vtkSmartPointer<vtkPolyData> dataTarget = vtkSmartPointer<vtkPolyData>::New();
    dataTarget = sourceTarget->GetPolyDataOutput();
//    dataTarget = sourceTarget->GetOutput();
    int pointCountTarget = dataTarget->GetNumberOfPoints();
    
    vtkSmartPointer<vtkPolyData> dataTemplate = vtkSmartPointer<vtkPolyData>::New();
    dataTemplate = sourceTemplate->GetPolyDataOutput();
//    dataTemplate = sourceTemplate->GetOutput();
    int pointCountTemplate = dataTemplate->GetNumberOfPoints();
    
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
    
    vtkSmartPointer<vtkPolyData> polydataTarget = vtkSmartPointer<vtkPolyData>::New();
    polydataTarget->SetPoints(pointsTarget);
    polydataTarget->SetPolys(dataTarget->GetPolys());
    
    vtkSmartPointer<vtkOctreePointLocator> octreeTarget = vtkSmartPointer<vtkOctreePointLocator>::New();
    octreeTarget->SetDataSet(polydataTarget);
    octreeTarget->BuildLocator();
    
    // Prepare to Find correspondences and the loop
    
    double maxDistance = 0;
    while (!maxDistance) {
        cout << "Please give the maximum distance between the corresponding points:" << endl;
        cin >> maxDistance;
    }
    
    double **finalDataTemplate;
    double distance;
    double overallDistance = 0;
    double minDistance = 1000000;
    int minOrientation = 0;
    int correspondences[pointCountTemplate];
    int currentCorrespondence;
    
    int discarded = 0;

    double epsilon = 1;
    
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
    
    vtkSmartPointer<vtkPolyData> polydataTemplate = vtkSmartPointer<vtkPolyData>::New();
    polydataTemplate->SetPoints(pointsTemplate);
    polydataTemplate->SetPolys(dataTemplate->GetPolys());
    
    vtkSmartPointer<vtkPolyData> transformedTemplate = vtkSmartPointer<vtkPolyData>::New();
    transformedTemplate->SetPoints(transformedPointsTemplate);
    transformedTemplate->SetPolys(dataTemplate->GetPolys());
    
    
//  A new vtkPoints object transformedPointsTemplate is done and filled with the final template points that are about to be displaced.
//  The scaling is done with the helper function void TransformScaleTranslate(vtkPoints points, vtkPolyData polydata, double value) with the
//  maxCoordinate variable and the three point sets:
//  template, target and transformed template. A visualization window is opened that displays the scaled meshes.
    
    // Scale to [-1,1] Cube
    
    polydataTarget = TransformScaleTranslate(pointsTarget, polydataTarget, 1/maxCoordinate);
    transformedTemplate = TransformScaleTranslate(transformedPointsTemplate, transformedTemplate, 1/maxCoordinate);
    polydataTemplate = TransformScaleTranslate(pointsTemplate, polydataTemplate, 1/maxCoordinate);
    
    pointsTarget = polydataTarget->GetPoints();
    transformedPointsTemplate = transformedTemplate->GetPoints();
    pointsTemplate = polydataTemplate->GetPoints();
    
    // Building the octree for the second time with the scaled mesh
    
    vtkSmartPointer<vtkPolyData> polydataResizedTargetOctree = vtkSmartPointer<vtkPolyData>::New();
    polydataResizedTargetOctree->SetPoints(pointsTarget);
    polydataResizedTargetOctree->SetPolys(polydataTarget->GetPolys());
    
    vtkSmartPointer<vtkOctreePointLocator> octreeResizedTarget = vtkSmartPointer<vtkOctreePointLocator>::New();
    octreeResizedTarget->SetDataSet(polydataResizedTargetOctree);
    octreeResizedTarget->BuildLocator();
    
    // Visualization
    
    vtkSmartPointer<vtkPolyDataMapper> polyDataMapperTarget = vtkSmartPointer<vtkPolyDataMapper>::New();
    polyDataMapperTarget->SetInput(polydataTarget);
    vtkSmartPointer<vtkPolyDataMapper> polyDataMapperTemplate = vtkSmartPointer<vtkPolyDataMapper>::New();
    polyDataMapperTemplate->SetInput(polydataTemplate);
    vtkSmartPointer<vtkPolyDataMapper> polyDataMapperTemplateTR = vtkSmartPointer<vtkPolyDataMapper>::New();
    polyDataMapperTemplateTR->SetInput(transformedTemplate);
    
    vtkSmartPointer<vtkProperty> targetProperty = vtkSmartPointer<vtkProperty>::New();
    targetProperty->SetOpacity(0.5);
    vtkSmartPointer<vtkProperty> templateProperty = vtkSmartPointer<vtkProperty>::New();
    templateProperty->SetColor(1, 0, 0);
    templateProperty->SetOpacity(0.2);
    vtkSmartPointer<vtkProperty> transformedTemplateProperty = vtkSmartPointer<vtkProperty>::New();
    transformedTemplateProperty->SetColor(0, 0, 1);
    
    vtkSmartPointer<vtkActor> actorTarget = vtkSmartPointer<vtkActor>::New();
    actorTarget->SetProperty(targetProperty);
    actorTarget->SetMapper(polyDataMapperTarget);
    actorTarget->SetPosition(0, 0, 0);
    
    vtkSmartPointer<vtkActor> actorTemplate = vtkSmartPointer<vtkActor>::New();
    actorTemplate->SetProperty(templateProperty);
    actorTemplate->SetMapper(polyDataMapperTemplate);
    actorTemplate->SetPosition(0, 0, 0);
    
    vtkSmartPointer<vtkActor> actorTemplateTr = vtkSmartPointer<vtkActor>::New();
    actorTemplateTr->SetProperty(transformedTemplateProperty);
    actorTemplateTr->SetMapper(polyDataMapperTemplateTR);
    actorTemplateTr->SetPosition(0, 0, 0);
    
    vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
    camera->SetPosition(1,1,1);
    camera->SetFocalPoint(0,0,0);
    
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actorTarget);
    //renderer->AddActor(actorTemplate);
    renderer->AddActor(actorTemplateTr);
    
    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
    renderer->SetBackground(1,1,1);
    
    vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
    window->SetSize(1200, 800);
    window->AddRenderer(renderer);
    
    window->Render();
    window->Start();
    
//  We enter the first loop and right after this the second.
//  The correspondences are set and the equation is build and solved with the help of SparseLU Eigen class.
//  (BiCGSTAB class implementation is commented)
//  The newly found transformations are applied and the norm is calculated and compared to the epsilon value.
    
    // Matrix X
    
    MatrixXd Xmatrix(4*pointCountTemplate, 3);
    Xmatrix.setZero();
    
    // The loop should start here
    
    double stiffness = 0;
    double minStiffness = 0;
    double stiffnessStep = 0;

    while (!minStiffness || !stiffness || !stiffnessStep) {
        cout << "Please give the maximum siffness value:" << endl;
        cin >> stiffness;
        cout << "And the minimum siffness value:" << endl;
        cin >> minStiffness;
        cout << "And the siffness step value:" << endl;
        cin >> stiffnessStep;
    }
    
    cout << "Overall points in the mesh: " << pointCountTemplate << endl << endl;
    cout << "Processing..." << endl;
    
    while (stiffness > minStiffness) {
    
        int loop = 1;
        while (loop) {
        
            for(int i = 0; i < pointCountTemplate; ++i)
            {
                double point[3];
                transformedPointsTemplate->GetPoint(i, point);
                currentCorrespondence = octreeResizedTarget->FindClosestPoint( point[0], point[1], point[2], distance );
                if (distance < maxDistance) {
                    correspondences[i] = currentCorrespondence;
                }
                else {
                    correspondences[i] = -1;
                    discarded++;
                }
            }
            
            cout << "Discarded points in this round of finding correspondences: " << discarded << endl;
            discarded = 0;
            
            // Equation
            
            // Extract edges
            
            vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
            extractEdges->SetInputConnection(sourceTemplate->GetOutputPort());
            extractEdges->Update();
            
            vtkCellArray* linesTemplate = extractEdges->GetOutput()->GetLines();
            int linesCountTemplate = linesTemplate->GetNumberOfCells();
            
            linesTemplate->InitTraversal();
            vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
            
            // Matrix A
            
            typedef Eigen::Triplet<double> T;
            vector<T> aTripletList;
            aTripletList.reserve(linesCountTemplate*8 + 4*pointCountTemplate);
            
            linesTemplate->GetNextCell(idList);
            for (int i = 0; i < linesCountTemplate; i++) {
                if (idList->GetId(0) < idList->GetId(1)) {
                    for (int j = 0; j < 4; j++) {
                        aTripletList.push_back(T(4*i + j, 4*idList->GetId(0) + j,-1*stiffness));
                        aTripletList.push_back(T(4*i + j, 4*idList->GetId(1) + j,1*stiffness));
                    }
                }
                else if (idList->GetId(1) < idList->GetId(0)) {
                    for (int j = 0; j < 4; j++) {
                        aTripletList.push_back(T(4*i + j, 4*idList->GetId(1) + j,-1*stiffness));
                        aTripletList.push_back(T(4*i + j, 4*idList->GetId(0) + j,1*stiffness));
                    }
                }
                linesTemplate->GetNextCell(idList);
            }
            
            for (int i = 0; i < pointCountTemplate; i++) {
                double point[3];
                int w = correspondences[i]!=-1 ? 1 : 0;
                transformedPointsTemplate->GetPoint(i, point);
                for (int j = 0; j < 3; j++) {
                    aTripletList.push_back(T(4*linesCountTemplate + i, i*4 + j, point[j]*w));
                }
                aTripletList.push_back(T(4*linesCountTemplate + i, i*4 + 3, 1*w));
            }

            SparseMatrix<double> Amatrix(4*linesCountTemplate + pointCountTemplate, 4*pointCountTemplate);
            Amatrix.setFromTriplets(aTripletList.begin(), aTripletList.end());
            
            SparseMatrix<double> AmatrixT(4*pointCountTemplate, 4*linesCountTemplate + pointCountTemplate);
            AmatrixT = Amatrix.transpose();
            SparseMatrix<double> AmatrixFinal(4*pointCountTemplate, 4*pointCountTemplate);
            AmatrixFinal = AmatrixT * Amatrix;
            
            // Matrix B
            
            MatrixXd Bmatrix(4*linesCountTemplate + pointCountTemplate, 3);
            Bmatrix.setZero();

            for (int i = 0; i < pointCountTemplate; i++) {
                double point[3];
                int w = correspondences[i]!=-1 ? 1 : 0;
                pointsTarget->GetPoint(correspondences[i], point);
                for (int j = 0; j < 3; j++) {
                    Bmatrix(4*linesCountTemplate + i, j) = point[j]*w;
                }
            }

            MatrixXd BmatrixFinal(4*pointCountTemplate, 3);
            BmatrixFinal.setZero();
            BmatrixFinal = AmatrixT * Bmatrix;
            
            // Temporal Matrix X
            
            MatrixXd TempXmatrix(4*pointCountTemplate, 3);
            TempXmatrix = Xmatrix;
            
            // Solver
            
            clock_t start;
            double duration;
            
            start = clock();
            
//            BiCGSTAB<SparseMatrix<double> > solver;
//            AmatrixFinal.makeCompressed();
//            solver.compute(AmatrixFinal);
//            Xmatrix = solver.solve(BmatrixFinal);
            
            SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver;
            solver.analyzePattern(AmatrixFinal);
            solver.factorize(AmatrixFinal);
            Xmatrix = solver.solve(BmatrixFinal);
            
            duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
            
            cout<< "Duration of this solver round: " << duration << " seconds" << endl;

            // Transformation
            
            for (int i = 0; i < pointCountTemplate; i++) {
                double point[3];
                double result[3];
                transformedPointsTemplate->GetPoint(i, point);
                
                for (int j = 0; j < 3; j++) {
                    result[j] = point[0] * Xmatrix(i*4 + 0,j) + point[1] * Xmatrix(i*4 + 1,j) + point[2] * Xmatrix(i*4 + 2,j) + 1 * Xmatrix(i*4 + 3,j);
                }
                transformedPointsTemplate->InsertPoint(i, result);
                
            }
            
//          The visualization window is updated.
//          The inner loop ends when the norm value is smaller than the epsilon value.
//          The outer loop ends after some iterations with different stiffness weights.
            
            double norm = (Xmatrix - TempXmatrix).norm();
            
            cout << "The norm value in this round (inner loop): "<< norm << endl << endl;
            
            //loop = 0;
            if (norm < epsilon || norm != norm) {
                loop = 0;
                if (norm != norm) {
                    cout << "Bad Matrix Error!" << endl;
                    errorInt = 1;
                }
            }
            
            window->Finalize();
            window->Render();
            window->Start();
        }
        cout << "The stiffness value in this round (outer loop): "<< stiffness << endl;
        cout << "======================================================== "<< endl << endl;
        stiffness -= stiffnessStep;
        if (errorInt) {
            break;
        }
        
    }
    
    window->Finalize();
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->SetSize(1200, 800);
    renWin->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renWin);
    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());

    interactor->Initialize();
    interactor->Start();
    
    freeData(finalDataTarget);
    freeData(finalDataTemplate);
    Xmatrix.resize(0,0);
    
    return 0;
}
