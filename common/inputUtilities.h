#ifndef __SurfaceRegistration__inputUtilities__
#define __SurfaceRegistration__inputUtilities__

#include <stdio.h>
#include <iostream>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>

vtkSmartPointer<vtkPolyData> vtkGetInput (std::string path);
vtkSmartPointer<vtkPolyData> vtpGetInput (std::string path);
vtkSmartPointer<vtkDataSetReader> vtkGetSource (std::string path);
vtkSmartPointer<vtkXMLPolyDataReader> vtpGetSource (std::string path);

#endif
