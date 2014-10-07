#ifndef __SurfaceRegistration__visualizationUtilities__
#define __SurfaceRegistration__visualizationUtilities__

#include <stdio.h>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTrackball.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>

int visualizeWindow(vtkSmartPointer<vtkPolyData> polydataTarget,
               vtkSmartPointer<vtkPolyData> polydataTemplate,
               vtkSmartPointer<vtkPolyData> transformedTemplate,
               vtkSmartPointer<vtkRenderWindow> window, vtkSmartPointer<vtkRenderer> renderer);

#endif
