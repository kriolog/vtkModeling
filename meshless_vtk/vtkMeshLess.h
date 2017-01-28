/***********************************************************

Filename: vtkMeshLess.h

Description:

This class is for meshless model. Basically it uses distance to find
neighbors within a meshless dataset. Then use Laplacian Surface
Editing or Continuum Elasticity to calculate displacements. The
Laplacian part is almost the same as class vtkLaplacianSurface. The
only difference is the definition of neighbors. So this class extends
vtkLaplacianSurface for convenience.

Using LSE with meshless data: proposed by Xiaoxu Wang in "LV Motion
and Strain Computation from tMRI based on Meshless Deformable Models"

Using Continuum Elasticity with meshless data: proposed by Matthias
Muller in "Point Based Animation of Elastic, Plastic and Melting
Objects"

Input: Polydata (no connectivity), control points' id and positions,
subdivision level (using spatial hashing to divide the space, in order
to find neighbors efficiently. subdivision level means "how many
pieces divided along the x axis"), LSE or CE, threshold (for stopping
criteria)

Output: Deformed polydata (no connectivity)

Usage: see build/bin/test_meshless or python/meshlesstest.py

Author: Shaoting Zhang, Aug. 2008

***********************************************************/

#ifndef __vtkMeshLess_h__
#define __vtkMeshLess_h__

#include <vtkPolyDataAlgorithm.h>
#include <vtkIdList.h>
#include <vector>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <lse_vtk/vtkLaplacianSurface.h>

class vtkLaplacianSurface;

class VTK_EXPORT vtkMeshLess : public vtkLaplacianSurface {
public:
    static vtkMeshLess* New();
    vtkTypeMacro(vtkMeshLess,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    //get parameters for meshless model
    vtkSetMacro(IsLaplacianVolume, int);	// Set to 1 means Enable Laplacian Editing, default is 0
    vtkSetMacro(IsContinuumElasticity, int);	// Set to 1 means Enable Continuum Elasticity, default is 0
    vtkSetMacro(Subdivision, int);              // How many pieces divided along x axis
    vtkSetMacro(PoissonRatio, double);          // Poisson's ratio (material property)
    vtkSetMacro(YoungModulus, double);          // Young's Modulous (material property)
    vtkSetMacro(DeltaT, double);                // dt for every step

protected:
    vtkMeshLess();
    ~vtkMeshLess();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    //data
    int IsLaplacianVolume;	// use Laplacian Surface (or Volume) Editing
    int IsContinuumElasticity;	// use Continuum Elasticity
    int Subdivision;            // How many pieces divided along x axis
    double YoungModulus;        // Material property
    double PoissonRatio;        // Material property
    double DeltaT;              // dt for every step

    //methods
    Matrix neighborMatrix(vtkPolyData* input);        // Get the neighbor list as a matrix
    void CE(vtkPolyData* input, vtkPolyData* output); // Using Continuum Elasticity

private:

    vtkMeshLess(const vtkMeshLess&);    // Not implemented.
    void operator=(const vtkMeshLess&); // Not implemented.
};

#endif

