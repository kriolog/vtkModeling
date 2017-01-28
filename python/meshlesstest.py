from vtk import *
from vtkExtend import *

r = vtkPolyDataReader()
r.SetFileName("meshless/sphere.vtk")

'''
p = vtkPointSource()
p.SetNumberOfPoints(801)
p.SetRadius(10)

w = vtkPolyDataWriter()
w.SetFileName("meshless/sphere.vtk")
w.SetInputConnection(p.GetOutputPort())
w.Write()
'''
test = vtkMeshLess()
test.SetInputConnection(r.GetOutputPort())
test.SetIsLaplacianSurface(1)
pts = vtkPoints()
pts.InsertNextPoint(1,0,0)
test.SetControlPoints(pts)
ids = vtkIdList()
ids.InsertNextId(0)
test.SetControlIds(ids)

mapper = vtkPolyDataMapper()
mapper.SetInputConnection(test.GetOutputPort())

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
