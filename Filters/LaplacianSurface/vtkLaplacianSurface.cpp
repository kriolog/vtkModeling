#include "vtkLaplacianSurface.h"

#include <math.h>
#include <iostream>
#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkType.h>
#include "VectorMath.h"

using namespace std;

vtkStandardNewMacro(vtkLaplacianSurface);

// Default constructor.
vtkLaplacianSurface::vtkLaplacianSurface() {
    this->ver = vtkPoints::New();
    this->verPre = vtkPoints::New();
    this->verModify = vtkPoints::New();
    this->polys = vtkCellArray::New();
    this->writer = vtkPolyDataWriter::New();
    this->IsEditing = 0;
    this->IsOptimization = 0;
    this->IsCotangent = 0;
    this->Threshold = 0.00001;
    this->Kernel = 1;
    this->DisplayFrames = 0;

    this->SetNumberOfInputPorts(2);
}

// Destrucor
vtkLaplacianSurface::~vtkLaplacianSurface() {
    this->ver->Delete();
    this->verPre->Delete();
    this->verModify->Delete();
    this->polys->Delete();
    this->writer->Delete();
}

vtkPointSet* vtkLaplacianSurface::GetControlPointsData()
{
  if (this->GetNumberOfInputConnections(1) < 1) {
    return NULL;
  }

  return vtkPointSet::SafeDownCast(this->GetInputDataObject(1, 0));
}

void vtkLaplacianSurface::RoutineCheck(vtkPolyData* input) {
    if (!ControlIds
        || !ControlPoints
        || this->ControlIds->GetNumberOfTuples() == 0
        || this->ControlPoints->GetNumberOfPoints() == 0)
    {
        cout << "Control Id or Control Point is empty" << endl;
        exit(0);
    }
    if (this->ControlIds->GetNumberOfTuples() !=
            this->ControlPoints->GetNumberOfPoints()) {
        cout << "The number of control points and control ids should be the same" << endl;
        exit(0);
    }
    if (input->GetNumberOfPoints() < this->ControlIds->GetNumberOfTuples()) {
        cout << "Control points are more then points in the data" << endl;
        exit(0);
    }
    if (input == NULL) {
        cout << "Input data is empty" << endl;
        exit(0);
    }
}

// VTK specific method: This method is called when the pipeline is calculated.
int vtkLaplacianSurface::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    // get the input and ouptut
    vtkPolyData *input = vtkPolyData::SafeDownCast(
                             inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
                              outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkPointSet* ControlPointSet = this->GetControlPointsData();
    this->ControlPoints = ControlPointSet->GetPoints();
    this->ControlIds
      = vtkIdTypeArray::SafeDownCast(this->GetInputArrayToProcess(0, ControlPointSet));

    // Do LSO for smoothing or do LSE for deformation
    if (this->IsOptimization == 1) {
        LSO(input, output);
    } else if (this->IsEditing == 1) {
        this->RoutineCheck(input);
        LSE(input, output);
    } else {
        cout << "Need to specify the algorithm: Laplacian Mesh Optimization or Laplacian Surface Editing" << endl;
        exit(0);
    }
    return 1;
}

void vtkLaplacianSurface::SetAlgorithm(int algorithm)
{
  if(algorithm == 0) {
    IsOptimization = 1;
    IsEditing = 0;
  } else {
    IsOptimization = 0;
    IsEditing = 1;
  }
}

// meaning of variable names, see paper "Laplacian Mesh Optimization", formula (10)
void vtkLaplacianSurface::LSO(vtkPolyData* input, vtkPolyData* output) {
    output->DeepCopy(input);
    this->ver->DeepCopy(input->GetPoints());
    Matrix L = LMatrix(input);
    int sizeIds = this->ControlIds->GetNumberOfTuples();
    int sizePts = L.RowNo();
    Matrix A(sizePts+sizeIds, L.ColNo());
    Matrix b(sizePts+sizeIds, 3);
    A.Null();
    b.Null();
    A.PartCopy(L);
    if (IsCotangent == 1) {
        Matrix delMatC = delCMatrix(input);
        b.PartCopy(delMatC);
    }
    for (int i=0; i<sizeIds; i++) {
        A(sizePts+i, this->ControlIds->GetValue(i)) = 1;
        double p[3];
        this->ver->GetPoint(this->ControlIds->GetValue(i), p);
        b.SetRow(sizePts+i, p);
    }

    Matrix AT = ~A;
    Matrix ATA = AT.SparseMul(A);
    Matrix ATb = AT.SparseMul(b);
    Matrix b1 = ATb.GetCol(0);
    Matrix b2 = ATb.GetCol(1);
    Matrix b3 = ATb.GetCol(2);
    Matrix vMat(sizePts,3);
    vMat.Null();
    vMat.SetCol(0,ATA.SSDSolve(b1));
    vMat.SetCol(1,ATA.SSDSolve(b2));
    vMat.SetCol(2,ATA.SSDSolve(b3));
    this->ver = matrix2Pts(vMat);
    output->SetPoints(this->ver);
}

// meaning of variable names, see paper "Laplacian Surface Editing" by Olga Sorkine, formula (6)
void vtkLaplacianSurface::LSE(vtkPolyData* input, vtkPolyData* output) {
    vector<Matrix> AMat;
    vector<Matrix> TMat;
    output->DeepCopy(input);
    this->verPre->DeepCopy(input->GetPoints());
    vtkIdType sizePts = this->verPre->GetNumberOfPoints();
    this->verModify->DeepCopy(input->GetPoints());
    // Move control points to specific positions (this->ControlPoints)
    for (int i=0; i<this->ControlIds->GetNumberOfTuples(); i++) {
        double pt[3];
        this->ControlPoints->GetPoint(i, pt);
        this->verModify->SetPoint(this->ControlIds->GetValue(i), pt[0], pt[1], pt[2]);
    }
    Matrix verPreMat = pts2Matrix(this->verPre);
    Matrix verModifyMat = pts2Matrix(this->verModify);
    Matrix neiMat = neighborMatrix(input);
    // Calculate matrix: inverse of ((A transpose) * A) * (A transpose)
    // Put all such matrices in a vector, will be used later
    for (int i=0; i<sizePts; i++) {
        Matrix A((neiMat(i,0)+1)*3, 7);
        A.Null();
        for (int j=0; j<neiMat(i,0); j++) {
            int ptId = neiMat(i, j+1);
            int AId = j*3;
            A(AId, 0)=verPreMat(ptId,0);
            A(AId, 2)=verPreMat(ptId,2);
            A(AId, 3)=-verPreMat(ptId,1);
            A(AId, 4)=1;
            A(AId+1, 0)=verPreMat(ptId,1);
            A(AId+1, 1)=-verPreMat(ptId,2);
            A(AId+1, 3)=verPreMat(ptId,0);
            A(AId+1, 5)=1;
            A(AId+2, 0)=verPreMat(ptId,2);
            A(AId+2, 1)=verPreMat(ptId,1);
            A(AId+2, 2)=-verPreMat(ptId,0);
            A(AId+2, 6)=1;
        }
        Matrix ATA = (~A)*A;
        Matrix ATAA = (ATA.Inv())*(~A);
        AMat.push_back(ATAA);
    }
    Matrix errMat = verPreMat-verModifyMat;
    double errPre = errMat.Norm() + this->Threshold;
    // The loop will stop if current error is bigger than previous error
    // or smaller than a threshold
    while (errPre > errMat.Norm() && abs(errMat.normI())>this->Threshold) {
        errPre = abs(errMat.Norm());
        // Calculate matrix b and T
        for (int i=0; i<sizePts; i++) {
            Matrix b((neiMat(i,0)+1)*3, 1);
            b.Null();
            for (int j=0; j<neiMat(i,0); j++) {
                int ptId = neiMat(i, j+1);
                int AId = j*3;
                b(AId,0)=verModifyMat(ptId,0);
                b(AId+1,0)=verModifyMat(ptId,1);
                b(AId+2,0)=verModifyMat(ptId,2);
            }
            Matrix sht = AMat.at(i)*b;
            Matrix T(4,4);
            T.Null();
            T(0,0)=sht(0,0);
            T(0,1)=-sht(3,0);
            T(0,2)=sht(2,0);
            T(0,3)=sht(4,0);
            T(1,0)=sht(3,0);
            T(1,1)=sht(0,0);
            T(1,2)=-sht(1,0);
            T(1,3)=sht(5,0);
            T(2,0)=-sht(2,0);
            T(2,1)=sht(1,0);
            T(2,2)=sht(0,0);
            T(2,3)=sht(6,0);
            T(3,3)=1;
            TMat.push_back(T);
        }
        // Apply T on every point
        for (int i=0; i<sizePts; i++) {
            Matrix verTmp(4,1);
            verTmp(0,0)=verPreMat(i,0);
            verTmp(1,0)=verPreMat(i,1);
            verTmp(2,0)=verPreMat(i,2);
            verTmp(3,0)=1;
            verTmp = TMat.at(i)*verTmp;
            verModifyMat(i,0)=verTmp(0,0);
            verModifyMat(i,1)=verTmp(1,0);
            verModifyMat(i,2)=verTmp(2,0);
        }
        // Write current data to disk for later displaying
        // If "DisplayFrames" is set to 1
        if (this->DisplayFrames > 0) {
            this->DisplayFrames++;
            vtkPolyData* tempGrids = vtkPolyData::New();
            tempGrids->DeepCopy(output);
            this->verModify = matrix2Pts(verModifyMat);
            tempGrids->SetPoints(this->verModify);
            this->writer->SetInputData(tempGrids);
            char fn[10];
            this->filename(this->DisplayFrames, fn);
            this->writer->SetFileName(fn);
            this->Modified();
            this->writer->Write();
            tempGrids->Delete();
        }
        // Get error matrix
        Matrix verModifyMatPre = verModifyMat;
        for (int i=0; i<this->ControlIds->GetNumberOfTuples(); i++) {
            double pt[3];
            this->ControlPoints->GetPoint(i, pt);
            verModifyMat(this->ControlIds->GetValue(i), 0)=pt[0];
            verModifyMat(this->ControlIds->GetValue(i), 1)=pt[1];
            verModifyMat(this->ControlIds->GetValue(i), 2)=pt[2];
        }
        errMat = verModifyMatPre - verModifyMat;
        TMat.clear();
    } // end of while
    this->verModify = matrix2Pts(verModifyMat);
    output->SetPoints(this->verModify);
    AMat.clear();
    TMat.clear();
}

Matrix vtkLaplacianSurface::neighborMatrix(vtkPolyData* input) {
    vtkIdType sizePts = input->GetNumberOfPoints();
    vtkIdType sizePly = input->GetNumberOfPolys();
    Matrix adjTemp(sizePts,sizePts);
    adjTemp.Null();
    Matrix adjMat(sizePts,sizePts);
    adjMat.Unit();
    this->polys = input->GetPolys();
    this->polys->InitTraversal();
    vtkIdType* cell;
    vtkIdType ncell;
    for (int i=0; i<sizePly; i++) {
        this->polys->GetNextCell(ncell, cell);
        adjTemp(cell[0], cell[1]) = 1;
        adjTemp(cell[1], cell[0]) = 1;
        adjTemp(cell[1], cell[2]) = 1;
        adjTemp(cell[2], cell[1]) = 1;
        adjTemp(cell[0], cell[2]) = 1;
        adjTemp(cell[2], cell[0]) = 1;
    }
    for (int i=0; i<this->Kernel; i++) {
        adjMat = adjMat * adjTemp;
    }
    Matrix neiMat(sizePts, adjMat.norm1()+2);
    neiMat.Null();
    for (int i=0; i<sizePts; i++) {
        int rowNum = 1;
        neiMat(i,rowNum) = i;
        for (int j=0; j< sizePts; j++) {
            if (adjMat(i,j) >= 1) {
                rowNum++;
                neiMat(i,rowNum) = j;
            }
        }
        neiMat(i,0) = rowNum;
    }
    return neiMat;
}

Matrix vtkLaplacianSurface::LMatrix(vtkPolyData* input) {
    vtkIdType sizePts = input->GetNumberOfPoints();
    vtkIdType sizePly = input->GetNumberOfPolys();
    Matrix adjMat(sizePts, sizePts);
    Matrix DMat(sizePts, sizePts);
    adjMat.Null();    // set all elements to zero
    DMat.Null();
    this->polys = input->GetPolys();
    this->polys->InitTraversal();
    vtkIdType* cell;
    vtkIdType ncell;

    for (int i=0; i<sizePly; i++) {
        polys->GetNextCell(ncell, cell);
        adjMat(cell[0], cell[1]) = 1;
        adjMat(cell[1], cell[0]) = 1;
        adjMat(cell[1], cell[2]) = 1;
        adjMat(cell[2], cell[1]) = 1;
        adjMat(cell[0], cell[2]) = 1;
        adjMat(cell[2], cell[0]) = 1;
        DMat(cell[0], cell[0]) += 1;
        DMat(cell[1], cell[1]) += 1;
        DMat(cell[2], cell[2]) += 1;
    }

    Matrix IMat(sizePts, sizePts);
    IMat.Unit();    // make identity matrix
    Matrix LMat = IMat - (DMat.DiagInv()).SparseMul(adjMat);
    return LMat;
}

Matrix vtkLaplacianSurface::delCMatrix(vtkPolyData* input) {
    vtkIdType sizePts = input->GetNumberOfPoints();
    vtkIdType sizePly = input->GetNumberOfPolys();
    vtkPoints* pts = input->GetPoints();
    Matrix vMat = pts2Matrix(pts);
    this->polys = input->GetPolys();
    this->polys->InitTraversal();
    vtkIdType* cell;
    vtkIdType ncell;
    Matrix weightC(sizePts, sizePts);
    Matrix delMatC(sizePts, 3);
    weightC.Null();
    delMatC.Null();
    //Get the delMatC according to "Generalized Barycentric Coordinates on
    //Irregular Polygons" by Mark Meyery and Haeyoung Leez
    for (int i=0; i<sizePly; i++) {
        double v21[3], v23[3], v12[3], v13[3], v31[3], v32[3], p1[3], p2[3], p3[3];
        this->polys->GetNextCell(ncell, cell);
        pts->GetPoint(cell[0], p1);
        pts->GetPoint(cell[1], p2);
        pts->GetPoint(cell[2], p3);
        VectorMath::pts2Vec(p1, p2, v21);
        VectorMath::pts2Vec(p2, p1, v12);
        VectorMath::pts2Vec(p3, p2, v23);
        VectorMath::pts2Vec(p2, p3, v32);
        VectorMath::pts2Vec(p3, p1, v13);
        VectorMath::pts2Vec(p1, p3, v31);
        double angle123 = acos(VectorMath::DotProduct(v21,v23)/
                               (VectorMath::Length(v21)*VectorMath::Length(v23)));
        double angle312 = acos(VectorMath::DotProduct(v12,v13)/
                               (VectorMath::Length(v12)*VectorMath::Length(v13)));
        double angle132 = acos(VectorMath::DotProduct(v31,v32)/
                               (VectorMath::Length(v31)*VectorMath::Length(v32)));
        double area12 = VectorMath::Area(v12);
        double area13 = VectorMath::Area(v13);
        double area23 = VectorMath::Area(v23);
        weightC(cell[0], cell[1]) = weightC(cell[0], cell[1]) + 1/(area12*tan(angle123));
        weightC(cell[1], cell[0]) = weightC(cell[1], cell[0]) + 1/(area12*tan(angle312));
        weightC(cell[0], cell[2]) = weightC(cell[0], cell[2]) + 1/(area13*tan(angle132));
        weightC(cell[2], cell[0]) = weightC(cell[2], cell[0]) + 1/(area13*tan(angle312));
        weightC(cell[1], cell[2]) = weightC(cell[1], cell[2]) + 1/(area23*tan(angle132));
        weightC(cell[2], cell[1]) = weightC(cell[2], cell[1]) + 1/(area23*tan(angle123));
    }
    Matrix mSum = weightC.sum();
    for (int i=0; i<sizePts; i++) {
        for (int j=0; j<sizePts; j++) {
            weightC(i,j) = weightC(i,j)/mSum(i,0);
        }
    }
    Matrix I(sizePts, sizePts);
    I.Unit();
    delMatC = (I - weightC).SparseMul(vMat);
    return delMatC;
}

vtkPoints* vtkLaplacianSurface::matrix2Pts(Matrix m) {
    vtkPoints* pts = vtkPoints::New();
    for (unsigned i=0; i<m.RowNo(); i++) {
        pts->InsertNextPoint(m(i,0), m(i,1), m(i,2));
    }
    return pts;
}

Matrix vtkLaplacianSurface::pts2Matrix(vtkPoints* v) {
    vtkIdType sizePts = v->GetNumberOfPoints();
    Matrix vMat(sizePts, 3);
    vMat.Null();
    for (int i=0; i<sizePts; i++) {
        double pt[3];
        v->GetPoint(i, pt);
        vMat(i,0)=pt[0];
        vMat(i,1)=pt[1];
        vMat(i,2)=pt[2];
    }
    return vMat;
}

void vtkLaplacianSurface::filename(int i, char* c) {
    sprintf(c, "%d", i);
}

// VTK specific method:
//      This method sets the input ports for this pipeline element.
int vtkLaplacianSurface::FillInputPortInformation(int port, vtkInformation* info) {
  if (port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  else
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");

  return 1;
}

// VTK specific method:
//      This method prints the contents of this object. The method is needed
//      to conform to the VTK tools.
void vtkLaplacianSurface::PrintSelf(ostream& os, vtkIndent indent) {
    vtkPolyDataAlgorithm::PrintSelf(os,indent);
}
