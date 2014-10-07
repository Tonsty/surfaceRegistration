#include "algorithm.h"

int mainAlgorithmLoop (double stiffness, double minStiffness, double stiffnessStep, double maxDistance,
                       vtkSmartPointer<vtkPoints> pointsTarget, vtkSmartPointer<vtkPoints> transformedPointsTemplate,
                       int pointCountTemplate, vtkSmartPointer<vtkOctreePointLocator> octreeResizedTarget,
                       vtkCellArray* linesTemplate,
                       vtkSmartPointer<vtkRenderWindow> window) {

    
    int correspondences[pointCountTemplate];
    int currentCorrespondence;
    
    double distance = 0;
    int discarded = 0;
    int errorInt = 0;
    double epsilon = 1;
    
    MatrixXd Xmatrix(4*pointCountTemplate, 3);
    Xmatrix.setZero();
    
    // The loop should start here
    
    std::cout << "Overall points in the mesh: " << pointCountTemplate << std::endl;
    std::cout << "Processing..." << std::endl;
    
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
            
            std::cout << "Discarded points in this round of finding correspondences: " << discarded << std::endl;
            discarded = 0;
            
            // Equation
            
            // Use the extract edges
            int linesCountTemplate = linesTemplate->GetNumberOfCells();
            
            linesTemplate->InitTraversal();
            vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
            
            // Matrix A
            
            typedef Eigen::Triplet<double> T;
            std::vector<T> aTripletList;
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
            
            std::cout<< "Duration of this solver round: " << duration << " seconds" << std::endl;
            
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
            
            std::cout << "The norm value in this round (inner loop): "<< norm << std::endl << std::endl;
            
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
        std::cout << "The stiffness value in this round (outer loop): "<< stiffness << std::endl;
        std::cout << "======================================================== "<< std::endl << std::endl;
        stiffness -= stiffnessStep;
        if (errorInt) {
            break;
        }
        
    }
    
    return 0;
}