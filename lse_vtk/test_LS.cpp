#include <iostream>
#include "vtkLaplacianSurface.h"
#include "vtkSphereSource.h"
#include "vtkDataSetMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkPolyDataWriter.h"

using namespace std;

int main(int argc, char* argv[]) {

    vtkSphereSource* sphere = vtkSphereSource::New();
    sphere->SetRadius(1.0);
    sphere->SetThetaResolution(30);
    sphere->SetPhiResolution(30);

    vtkPoints *pts = vtkPoints::New();
    vtkIdList *ids = vtkIdList::New();
    pts->InsertNextPoint(0,0,1.3);
    ids->InsertNextId(0);
    ids->InsertNextId(1);

    vtkLaplacianSurface *test = vtkLaplacianSurface::New();
    test->SetInput(sphere->GetOutput());
    test->SetControlPoints(pts);
    test->SetControlIds(ids);
    test->SetIsEditing(1);
    test->SetKernel(1);
    /*
    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetInput(sphere->GetOutput());
    w->SetFileName("test.vtk");
    w->Write();
    */
    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputConnection(test->GetOutputPort());

    vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);

    vtkRenderer *ren = vtkRenderer::New();
    ren->AddActor(actor);
    ren->SetBackground(0.1,0.2,0.4);

    vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer( ren );
    renWin->SetSize( 300, 300 );

    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

    vtkInteractorStyleTrackballCamera *style =
        vtkInteractorStyleTrackballCamera::New();
    iren->SetInteractorStyle(style);

    iren->Initialize();
    iren->Start();

    sphere->Delete();
    mapper->Delete();
    actor->Delete();
    ren->Delete();
    renWin->Delete();
    iren->Delete();
    style->Delete();
}
