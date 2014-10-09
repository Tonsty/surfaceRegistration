#ifndef __SurfaceRegistration__preProcessing__
#define __SurfaceRegistration__preProcessing__

#include <stdio.h>
#include <iostream>
#include "common/PCAStatistics.h"

#include <vtkOctreePointLocator.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

vtkSmartPointer<vtkPolyData> TransformScaleTranslate(vtkSmartPointer<vtkPolyData> polydata, double value);

vtkSmartPointer<vtkOctreePointLocator> buildOctree(vtkSmartPointer<vtkPolyData> newdata, vtkSmartPointer<vtkPoints> points,
                                                   vtkSmartPointer<vtkPolyData> olddata);

double applyPCA(int pointCountTarget, int pointCountTemplate,
                                vtkSmartPointer<vtkPolyData> dataTarget, vtkSmartPointer<vtkPolyData> dataTemplate,
                                vtkSmartPointer<vtkPolyData> polydataTarget, vtkSmartPointer<vtkPolyData> polydataTemplate, vtkSmartPointer<vtkPolyData> transformedTemplate,
                                int maxDistance);

#endif
