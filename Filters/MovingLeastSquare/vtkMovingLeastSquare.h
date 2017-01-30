/**************************************************************

Filename: vtkMovingLeastSquare.h

Description:

Implemented paper "Image Deformation Using Moving Least Squares"
Formula (5)-(8), Affine/Similarity/Rigid transformation

Essential inputs: control points' id, control points' coordinates,
which transformation (A/S/R)

Optional inputs: alpha value for computing weight (default is 1.5),
dimension (default is 2D)

Usage: see python/mlstest.py

Author: Shaoting Zhang

Notes: can only deal with 2D case now

**************************************************************/

#ifndef __vtkMovingLeastSquare_h__
#define __vtkMovingLeastSquare_h__

#include <vtkPolyDataAlgorithm.h>
#include <vtkIdTypeArray.h>

// For Matrix TCL, c++ lib for matrix operations
// http://www.techsoftpl.com/matrix/download.htm
#include "matrix.h"
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

class VTK_EXPORT vtkMovingLeastSquare : public vtkPolyDataAlgorithm {
public:
    static vtkMovingLeastSquare* New();
    vtkTypeMacro(vtkMovingLeastSquare,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Get parameters
    vtkSetMacro(IsAffine, int);
    vtkSetMacro(IsSimilarity, int);
    vtkSetMacro(IsRigid, int);
    vtkSetMacro(Dimension, int);
    vtkSetMacro(Alpha, double);

    vtkPointSet* GetControlPointsData();

protected:
    vtkMovingLeastSquare();
    ~vtkMovingLeastSquare();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int FillInputPortInformation(int, vtkInformation*);

    // data
    vtkIdTypeArray* Ids;      // control points' id
    vtkPoints* Q;        // control points' coordinates (after deformation)
    vtkPoints* P;        // control points' coordinates (before deformation)
    vtkPoints* pts;      // All points
    Matrix* PMat;        // Matrix representation for P
    Matrix* QMat;        // Matrix representation for Q
    Matrix* ptMat;       // Matrix representation for pts
    int IsAffine;        // Using Affine transformation if 1
    int IsSimilarity;    // Using Similarity transformation if 1
    int IsRigid;         // Using Rigid transformation if 1
    int Dimension;       // 2D or 3D, default value is 2
    double Alpha;        // Alpha value for weight. larger then 1 is a better choice

    // method
    Matrix affine(Matrix v, Matrix* p, Matrix* q);            // affine transformation
    Matrix similarity(Matrix v, Matrix* p, Matrix* q);        // similarity transformation
    Matrix simOper(Matrix m);                                 // a Matrix operator for similarity transformation
    Matrix rigid(Matrix v, Matrix* p, Matrix* q);             // rigid transformation
    double weight(Matrix pi, Matrix v);                       // Get weight (for moving least square)
    Matrix* pts2Matrix(vtkPoints* v);                         // Convert vtkPoints to Matrix
    vtkPoints* matrix2Pts(Matrix* m);                         // Convert Matrix to vtkPoints

private:
    vtkMovingLeastSquare(const vtkMovingLeastSquare&);  // Not implemented.
    void operator=(const vtkMovingLeastSquare&);        // Not implemented.
};

#endif

