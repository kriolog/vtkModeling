from vtk import *
from vtkExtend import *

##cube = vtkCubeSource()
##tri = vtkTriangleFilter()
##tri.SetInputConnection(cube.GetOutputPort())

sphere = vtkPolyDataReader()
sphere.SetFileName("lse/sphere.vtk")

#ptsRd = vtkPolyDataReader()
#ptsRd.SetFileName('lse/cp.vtk')
#ptsRd.Update()

pts = vtkPoints()
pts.InsertNextPoint(0,0,1.3)
#pts = ptsRd.GetOutput().GetPoints()

ids = vtkIdList()
ids.InsertNextId(0)
#f=open('lse/cpid.vtk', 'r')
#for line in f:
#    ids.InsertNextId(int(line))

test = vtkLaplacianSurface()
test.SetInput(sphere.GetOutput())
test.SetControlPoints(pts)
test.SetControlIds(ids)
test.SetIsEditing(1)

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

##f.close()

iren.Start()
