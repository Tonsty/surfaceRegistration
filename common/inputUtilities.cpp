#include "inputUtilities.h"

vtkSmartPointer<vtkPolyData> vtkGetInput (std::string path) {
    const char* nameInput;
    
    nameInput = path.c_str();
    
    vtkSmartPointer<vtkDataSetReader> source = vtkSmartPointer<vtkDataSetReader>::New();
    
    source->SetFileName(nameInput);
    source->Update();
    
    return source->GetPolyDataOutput();
}

vtkSmartPointer<vtkPolyData> vtpGetInput (std::string path) {
    const char* nameInput;
    
    nameInput = path.c_str();
    
    vtkSmartPointer<vtkXMLPolyDataReader> source = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    
    source->SetFileName(nameInput);
    source->Update();
    
    return source->GetOutput();
    
}

vtkSmartPointer<vtkDataSetReader> vtkGetSource (std::string path) {
    const char* nameInput;
    
    nameInput = path.c_str();
    
    vtkSmartPointer<vtkDataSetReader> source = vtkSmartPointer<vtkDataSetReader>::New();
    
    source->SetFileName(nameInput);
    source->Update();
    
    return source;
}

vtkSmartPointer<vtkXMLPolyDataReader> vtpGetSource (std::string path) {
    const char* nameInput;
    
    nameInput = path.c_str();
    
    vtkSmartPointer<vtkXMLPolyDataReader> source = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    
    source->SetFileName(nameInput);
    source->Update();
    
    return source;
}