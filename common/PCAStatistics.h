#ifndef __SurfaceRegistration__PCAStatistics__
#define __SurfaceRegistration__PCAStatistics__

#include <stdio.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <vtkDoubleArray.h>
#include <vtkPCAStatistics.h>
#include <vtkTable.h>

double ** PCAStatistics(vtkSmartPointer<vtkPolyData> data, int pointCount, int orientation);

#endif
