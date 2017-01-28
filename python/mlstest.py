from vtk import *
from vtkExtend import *

mesh = vtkPolyDataReader()
mesh.SetFileName("mls/mesh.vtk")
mesh.Update()

pts = vtkPoints()
allpts = mesh.GetOutput().GetPoints()

ids = vtkIdList()
for i in range(0,10):
    ids.InsertNextId(i)
    pt = allpts.GetPoint(i)
    pts.InsertNextPoint(pt[0]-1,pt[1]-1,pt[2])
for i in range(90,100):
    ids.InsertNextId(i)
    pt = allpts.GetPoint(i)
    pts.InsertNextPoint(pt[0]+1,pt[1]+1,pt[2])
    
test = vtkMovingLeastSquare()
test.SetInput(mesh.GetOutput())
test.SetDimension(2)
test.SetIsRigid(1)
test.SetIds(ids)
test.SetQ(pts)

mapper = vtkPolyDataMapper()
mapper.SetInputConnection(test.GetOutputPort())

w = vtkPolyDataWriter()
w.SetInput(test.GetOutput())
w.SetFileName("mls/newMesh.vtk")
w.Write()

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
