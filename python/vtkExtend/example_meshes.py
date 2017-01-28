import vtk
import math
from __wrapped_classes import vtkGridVolumeMeshFilter


################################################################################
## generator for a meshed kidney
################################################################################
def kidney_surface(size = 1):
    kidney = vtk.vtkPolyDataReader()
    kidney.SetFileName("../../data/kidney_20.pdw")
    triangulated = vtk.vtkTriangleFilter()
    triangulated.SetInput(kidney.GetOutput())
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetTolerance(0.0001)
    cleaner.SetInput(triangulated.GetOutput())
    return cleaner


def kidney_volume(size = 1):
    kidney = vtk.vtkPolyDataReader()
    kidney.SetFileName("../../data/kidney_20.pdw")
    triangulated = vtk.vtkTriangleFilter()
    triangulated.SetInput(kidney.GetOutput())
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetTolerance(0.0001)
    cleaner.SetInput(triangulated.GetOutput())
    grid = vtkGridVolumeMeshFilter()
    grid.SetInputConnection(cleaner.GetOutputPort())
    grid.SetInteriorPointSpacing(10)
    return grid

################################################################################
## generator for a meshed cube
################################################################################
def cube(resolution = 7, size = 1):
    return cuboid(resolution, [size]*3)
    
def cuboid(resolution = 7, size = [1,1,1]):
    points = vtk.vtkPoints()
    
    dx = float(size[0])/float(resolution)
    dy = float(size[1])/float(resolution)
    dz = float(size[2])/float(resolution)
    
    xh = float(size[0])/2.0
    yh = float(size[1])/2.0
    zh = float(size[2])/2.0
    
    for x in range(0, resolution + 1):
        for y in range(0, resolution + 1):
            for z in range(0, resolution + 1):
                points.InsertNextPoint(x * dx - xh, y * dy - yh, z * dz - zh)
    
    grid = vtk.vtkPolyData()
    grid.SetPoints(points)
    
    mesh = vtk.vtkDelaunay3D()
    mesh.SetInput(grid)
    
    return mesh


################################################################################
## generator for a single tetrahedron
################################################################################
def tetra(size = 1):
    size = size / 2.0
    
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, -size, -size)
    points.InsertNextPoint(-size, -size, size)
    points.InsertNextPoint(size, -size, size)
    points.InsertNextPoint(0, size, 0)

    
    tetra = vtk.vtkTetra()
    tetra.GetPointIds().SetId(0,0)
    tetra.GetPointIds().SetId(1,1)
    tetra.GetPointIds().SetId(2,2)
    tetra.GetPointIds().SetId(3,3)


    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())
    
    mesh = vtk.vtkDelaunay3D()# vtk.vtkPassThroughFilter()
    mesh.SetInput(grid)
    
    return mesh


################################################################################
## generator for a cone
################################################################################
# todo



################################################################################
## generator for a meshed cylinder
################################################################################
def cylinder(resolution = 4, size = 1, ratio = 0.3):
    if resolution < 2:
        resolution = 2
    
    resolution *= 2
    
    grid = vtk.vtkAppendPolyData()
    
    # how to make facy cake:
    # add layers of pies... :-)
    height = float(size) / (resolution * 2 - 1)
    for l in range(resolution):
        # insert centers
        pts = vtk.vtkPoints()
        pts.InsertNextPoint(0, height * float((l - (resolution * 2 - 1) / 2) * 2 + 1.5), 0)
        pts.InsertNextPoint(0, height * float((l - (resolution * 2 - 1) / 2) * 2 + .5), 0)
        
        pol = vtk.vtkPolyData()
        pol.SetPoints(pts)
        
        grid.AddInput(pol)
        
        sides = 6
        for r in range(resolution / 2 - 1):
            pie = vtk.vtkCylinderSource()
            pie.SetHeight(height)
            pie.SetRadius(float(r + 1) * float(size) / 2.0 * float(ratio) / float(resolution/2 - 1))
            pie.SetCenter(0, height * float((l - (resolution * 2 - 1) / 2) * 2 + 1), 0)
            pie.SetResolution(sides)
            
            grid.AddInput(pie.GetOutput())
            
            sides *= 2
            
    
    mesh = vtk.vtkDelaunay3D()
    mesh.SetInput(grid.GetOutput())
    
    return mesh


################################################################################
## generator for a meshed sphere
################################################################################
def sphere(resolution = 3, size = 1):
    grid = vtk.vtkAppendPolyData()
    
    # add the center to the mesh
    center = vtk.vtkPolyData()
    center.SetPoints(vtk.vtkPoints())
    center.GetPoints().InsertNextPoint(0, 0, 0)
    
    grid.AddInput(center)
    
    # add spheres of increasing complexity
    subdivisions = 0
    radius = 0
    for i in range(resolution - 1):
        radius += float(size) / float(resolution - 1)
        
        sphere = nice_sphere(radius / 2.0, subdivisions)
        grid.AddInput(sphere)
        
        subdivisions += 1
    
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInput(grid.GetOutput())
    
    mesh = vtk.vtkDelaunay3D()
    mesh.SetInput(cleaner.GetOutput())
    
    return mesh


    
    
    
    
    
    
    
    
    
    
    
    
    
    
# some helper functions
###################################
def unit(v):
    l = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    return (v[0] / l, v[1] / l, v[2] / l)


def subdivide(a, b, c, depth):
    if depth == 0:
        return []
    d = unit((b[0] + c[0], b[1] + c[1], b[2] + c[2]))
    e = unit((a[0] + c[0], a[1] + c[1], a[2] + c[2]))
    f = unit((a[0] + b[0], a[1] + b[1], a[2] + b[2]))
    if depth > 1:
        return [d, e, f] + subdivide(a, f, e, depth - 1) \
                         + subdivide(d, e, f, depth - 1) \
                         + subdivide(f, b, d, depth - 1) \
                         + subdivide(e, d, c, depth - 1)
    else:
        return [d, e, f]




def nice_sphere(radius, subdivisions, center = (0, 0, 0)):
    platonic_solid = vtk.vtkPlatonicSolidSource()
    platonic_solid.SetSolidTypeToIcosahedron()
    platonic_solid = platonic_solid.GetOutput()
    platonic_solid.Update()
    
    points = []
    for i in range(platonic_solid.GetNumberOfCells()):
        cell = platonic_solid.GetCell(i)
        a = platonic_solid.GetPoint(cell.GetPointId(0))
        b = platonic_solid.GetPoint(cell.GetPointId(1))
        c = platonic_solid.GetPoint(cell.GetPointId(2))
        
        points += [a, b, c] \
               + subdivide(a, b, c, subdivisions)
        
    
    vtkPoints = vtk.vtkPoints()
    for point in points:
        vtkPoints.InsertNextPoint(center[0] + point[0] * radius,
                                  center[1] + point[1] * radius, 
                                  center[2] + point[2] * radius)
    poly = vtk.vtkPolyData()
    poly.SetPoints(vtkPoints)
    cells = vtk.vtkCellArray()
    for i in range(len(points)):
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, i)
        cells.InsertNextCell(vertex)
    poly.SetVerts(cells)
    return poly





