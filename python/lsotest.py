from vtk import *
from vtkExtend import *
from math import *

r = vtkPolyDataReader()

r.SetFileName("lso/LV_boundary.vtk")
r.Update()
ply = vtkPolyData()
ptsply = vtkPoints()
ply = r.GetOutput()
ptsply = ply.GetPoints()

# LSO for smooth
ids = vtkIdList()
for i in range(0,40):
    ids.InsertNextId(i*70)

lso = vtkLaplacianSurface()
lso.SetInput(ply)
lso.SetControlIds(ids)
lso.SetIsOptimization(1)
lso.Update()
w = vtkPolyDataWriter()
w.SetInput(lso.GetOutput())
w.SetFileName("lso/LV_lso.vtk")
w.Write()

mapper = vtkPolyDataMapper()
mapper.SetInputConnection(lso.GetOutputPort())

actor = vtkActor()
actor.SetMapper(mapper)

ren = vtkRenderer()
ren.AddActor(actor)

renWin = vtkRenderWindow()
renWin.AddRenderer(ren)

iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.SetBackground(0, 0, 0)
renWin.SetSize(600, 250)
iren.Initialize()
renWin.Render()
iren.Start()
