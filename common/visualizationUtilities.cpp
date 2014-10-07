#include "visualizationUtilities.h"

int visualizeWindow(vtkSmartPointer<vtkPolyData> polydataTarget,
               vtkSmartPointer<vtkPolyData> polydataTemplate,
               vtkSmartPointer<vtkPolyData> transformedTemplate,
               vtkSmartPointer<vtkRenderWindow> window, vtkSmartPointer<vtkRenderer> renderer) {

    vtkSmartPointer<vtkPolyDataMapper> polyDataMapperTarget = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkPolyDataMapper> polyDataMapperTemplate = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkPolyDataMapper> polyDataMapperTemplateTR = vtkSmartPointer<vtkPolyDataMapper>::New();

#if VTK_MAJOR_VERSION <= 5
    polyDataMapperTarget->SetInput(polydataTarget);
    polyDataMapperTemplate->SetInput(polydataTemplate);
    polyDataMapperTemplateTR->SetInput(transformedTemplate);
#else
    polyDataMapperTarget->SetInputData(polydataTarget);
    polyDataMapperTemplate->SetInputData(polydataTemplate);
    polyDataMapperTemplateTR->SetInputData(transformedTemplate);
#endif
    
    vtkSmartPointer<vtkProperty> targetProperty = vtkSmartPointer<vtkProperty>::New();
    targetProperty->SetOpacity(0.5);
    targetProperty->SetColor(1, 0, 0);
    vtkSmartPointer<vtkProperty> templateProperty = vtkSmartPointer<vtkProperty>::New();
    templateProperty->SetColor(0, 1, 0);
    templateProperty->SetOpacity(0.1);
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
    
    renderer->AddActor(actorTarget);
    //renderer->AddActor(actorTemplate);
    renderer->AddActor(actorTemplateTr);
    
    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();
    renderer->ResetCameraClippingRange();
    renderer->SetBackground(1,1,1);
    
    window->SetSize(1200, 800);
    window->AddRenderer(renderer);
    
    window->Render();
    window->Start();
    
    return 0;
}