#include <iostream>
#include "vtkMeshLess.h"
#include "vtkPointSource.h"
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
#include "vtkPolyDataReader.h"

using namespace std;

int main(int argc, char* argv[]) {
    vtkPolyDataReader *r = vtkPolyDataReader::New();
    r->SetFileName("sphere.vtk");
    r->Update();

    vtkPoints *ptssph = (r->GetOutput())->GetPoints();
    ptssph->Modified();

    vtkPoints *pts = vtkPoints::New();
    vtkIdList *ids = vtkIdList::New();
    for (int i=0; i<100; i++) {
        double pt[3];
        ptssph->GetPoint(i, pt);
        pt[0]=pt[0]+20;
        pt[1]=pt[1]+20;
        pt[2]=pt[2]+20;
        pts->InsertNextPoint(pt);
        ids->InsertNextId(i);
    }

    vtkMeshLess *test = vtkMeshLess::New();
    test->SetInput(r->GetOutput());
    test->SetControlPoints(pts);
    test->SetControlIds(ids);
    test->SetIsLaplacianVolume(1);

    vtkPolyDataWriter *w = vtkPolyDataWriter::New();
    w->SetInputConnection(test->GetOutputPort());
    w->SetFileName("test.vtk");
    w->Write();

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

    r->Delete();
    mapper->Delete();
    actor->Delete();
    ren->Delete();
    renWin->Delete();
    iren->Delete();
    style->Delete();
}
