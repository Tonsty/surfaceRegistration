#ifndef __SurfaceRegistration__algorithm__
#define __SurfaceRegistration__algorithm__

#include <stdio.h>

#include <vtkSmartPointer.h>
#include <vtkExtractEdges.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkRenderWindow.h>
#include <vtkOctreePointLocator.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/SparseExtra>

using namespace Eigen;

int mainAlgorithmLoop (double stiffness, double minStiffness, double stiffnessStep, double maxDistance,
                       vtkSmartPointer<vtkPoints> pointsTarget, vtkSmartPointer<vtkPoints> transformedPointsTemplate,
                       int pointCountTemplate, vtkSmartPointer<vtkOctreePointLocator> octreeResizedTarget,
                       vtkCellArray* linesTemplate,
                       vtkSmartPointer<vtkRenderWindow> window);

#endif
