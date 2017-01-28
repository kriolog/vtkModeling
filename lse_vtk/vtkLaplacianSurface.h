/***********************************************************

Filename: vtkLaplacianSurface.h

Description:

This class is for Laplacian Surface Editing/Optimization algorithms.

For LSE the inputs are control points id and positions. The class will
compute the deformation of rest points. Proposed by Olga Sorkine in
"Laplacian Surface Editing"

For LSO the inputs are just anchor points id. The class will optimize
rest triangles. Proposed by Andrew Nealen in "Laplacian Mesh
Optimization"

Usage: see python/lsetest.py and python/regis.py

Author: Shaoting Zhang

***********************************************************/

#ifndef __vtkLaplacianSurface_h__
#define __vtkLaplacianSurface_h__

#include <vtkPolyDataAlgorithm.h>
#include <vtkIdList.h>
#include <vector>
#include <vtkPolyDataWriter.h>

// For Matrix TCL, c++ lib for matrix operations
// http://www.techsoftpl.com/matrix/download.htm
#include "include/matrix.h"
#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
typedef matrix<int> iMatrix;
#else
typedef matrix Matrix;
#endif
// end of Matrix TCL

class vtkPolyData;

class VTK_EXPORT vtkLaplacianSurface : public vtkPolyDataAlgorithm {
public:
    static vtkLaplacianSurface* New();
    vtkTypeMacro(vtkLaplacianSurface,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    //get parameters for lse
    vtkSetObjectMacro(ControlIds, vtkIdList);	// Input, control points' id
    vtkSetObjectMacro(ControlPoints, vtkPoints);// Input, control points' coordinates
    vtkSetMacro(IsEditing, int);		// Set to 1 means Enable Laplacian Editing, default is 0
    vtkSetMacro(IsOptimization, int);		// Set to 1 means Enable Laplacian Optimization, default is 0
    vtkSetMacro(IsCotangent, int);		// Set to 1 means Enable cotangent weight, default is 0
    vtkSetMacro(Threshold, double);             // Stopping criterion
    vtkSetMacro(Kernel, int);
    vtkSetMacro(DisplayFrames, int);            // Set to 1 means write each frame to disk (for displaying later)

protected:
    vtkLaplacianSurface();
    ~vtkLaplacianSurface();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int FillInputPortInformation(int, vtkInformation*);

    //data
    vtkIdList* ControlIds;   // control points' id
    vtkPoints* ControlPoints;// control points' coordinates
    int IsEditing;	     // use LSE
    int IsOptimization;	     // use LSO
    int IsCotangent;	     // use cotangent weight (for LSO)
    double Threshold;	     // stopping criterion, used in LSE
    int Kernel;              // used in "neighborMatrix" as the distance
    int DisplayFrames;       // Whether or not save each frame
    vtkPoints* ver;	     // vertices of input, used in LSO
    vtkPoints* verPre;	     // vertices of input in previous step, used in LSE
    vtkPoints* verModify;    // vertices of input in this step, used in LSE
    vtkCellArray* polys;     // polygons of input
    vtkPolyDataWriter* writer;

    //methods
    void LSO(vtkPolyData* input, vtkPolyData* output);		// Laplacian Surface Optimization
    void LSE(vtkPolyData* input, vtkPolyData* output);		// Laplacian Surface Editing
    virtual Matrix neighborMatrix(vtkPolyData* input);		// Get the neighbor list as a matrix
    Matrix LMatrix(vtkPolyData* input);				// Get the L matrix (see Laplacian Surface Editing): L = I - Inv(D)*adjMat
    Matrix delCMatrix(vtkPolyData* input);			// Get the laplacian/delta matrix with cotangent weights: = (I - weightC)*v
    Matrix pts2Matrix(vtkPoints* v);				// Convert vtkPoints to Matrix
    vtkPoints* matrix2Pts(Matrix m);				// Convert Matrix to vtkPoints
    void RoutineCheck(vtkPolyData* input);

    void filename(int i, char* c);    // input is the current frame "i", output is the corresponding filename "c"

private:

    vtkLaplacianSurface(const vtkLaplacianSurface&);    // Not implemented.
    void operator=(const vtkLaplacianSurface&);         // Not implemented.
};

#endif

