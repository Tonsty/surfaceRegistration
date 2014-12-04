#include "inputUtilities.h"

BGLUndirectedGraph vtkGetInput (std::string path) {
    const char* nameInput;
    
    nameInput = path.c_str();
    
    vtkSmartPointer<vtkDataSetReader> source = vtkSmartPointer<vtkDataSetReader>::New();
    
    source->SetFileName(nameInput);
    source->Update();
    
    const std::size_t V = source->GetPolyDataOutput()->GetNumberOfPoints();
    
    vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
    extractEdges->SetInputConnection(source->GetOutputPort());
    extractEdges->Update();
    
    const std::size_t E = extractEdges->GetOutput()->GetNumberOfCells();
    
    BGLUndirectedGraph g(V);
    
    for(int i = 0; i < E; i++)
    {
        vtkSmartPointer<vtkLine> line = vtkLine::SafeDownCast(extractEdges->GetOutput()->GetCell(i));
        add_edge(line->GetPointIds()->GetId(0), line->GetPointIds()->GetId(1), g);
    }
    
    return g;
}

BGLUndirectedGraph vtpGetInput (std::string path) {
    const char* nameInput;
    
    nameInput = path.c_str();
    
    vtkSmartPointer<vtkXMLPolyDataReader> source = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    
    source->SetFileName(nameInput);
    source->Update();
    
    const std::size_t V = source->GetOutput()->GetNumberOfPoints();
    
    vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
    extractEdges->SetInputConnection(source->GetOutputPort());
    extractEdges->Update();
    
    const std::size_t E = extractEdges->GetOutput()->GetNumberOfCells();
    
    BGLUndirectedGraph g(V);
    
    for(int i = 0; i < E; i++)
    {
        vtkSmartPointer<vtkLine> line = vtkLine::SafeDownCast(extractEdges->GetOutput()->GetCell(i));
        add_edge(line->GetPointIds()->GetId(0), line->GetPointIds()->GetId(1), g);
    }
    
    return g;
    
}