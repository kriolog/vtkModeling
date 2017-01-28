from vtk import *
from vtkExtend import *

mesh = vtkUnstructuredGridReader()
mesh.SetFileName("mssAnimation/cube.ugrid")
mesh.Update()

pts = vtkPoints()
ids = vtkIdList()
allpts = mesh.GetOutput().GetPoints()
for i in range(0,6):
    ids.InsertNextId(i)
    pt = allpts.GetPoint(i)
    pts.InsertNextPoint(pt[0]-.1,pt[1]-.1,pt[2]-.1)
for i in range(210,216):
    ids.InsertNextId(i)
    pt = allpts.GetPoint(i)
    pts.InsertNextPoint(pt[0]+.1,pt[1]+.1,pt[2]+.1)
    
test = vtkMassSpring()
test.SetInput(mesh.GetOutput())
test.SetControlIds(ids)
test.SetControlPoints(pts)
test.SetDistanceCoefficient(20)
test.SetDampingCoefficient(5)
test.SetMass(9)
test.SetDeltaT(0.001)
test.SetSteps(4000)
test.SetDisplayFrames(100)

mapper = vtkDataSetMapper()
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
