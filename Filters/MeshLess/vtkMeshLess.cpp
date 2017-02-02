#include "vtkMeshLess.h"

#define _USE_MATH_DEFINES
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

vtkStandardNewMacro(vtkMeshLess);


// Default constructor.
vtkMeshLess::vtkMeshLess() {
    this->IsLaplacianVolume = 0;
    this->IsContinuumElasticity = 0;
    this->Subdivision = 10;
    this->DeltaT = 0.01;
    this->PoissonRatio = 0.5;
    this->YoungModulus = 0.05;

    this->SetNumberOfInputPorts(Superclass::GetNumberOfInputPorts());
}

// Destrucor
vtkMeshLess::~vtkMeshLess() {
}

void vtkMeshLess::SetAlgorithm(int algorithm)
{
  if(algorithm == 0) {
    IsLaplacianVolume = 1;
    IsContinuumElasticity = 0;
  } else {
    IsLaplacianVolume = 0;
    IsContinuumElasticity = 1;
  }
  this->Modified();
}

// VTK specific method: This method is called when the pipeline is calculated.
int vtkMeshLess::RequestData(
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

    // Validate the input data
    this->RoutineCheck(input);
    if (this->Subdivision < 2) {
        cout << "Subdivision should be bigger than or equal to 2" << endl;
        exit(0);
    }
    // Using Laplacian Volume Editing or Continuum Elasticity
    if (this->IsLaplacianVolume == 1) {
        LSE(input, output);
    } else if (this->IsContinuumElasticity == 1) {
        CE(input, output);
    } else {
        cout << "Need to specify IsLaplacianVolume or IsContinuumElasticity" << endl;
        exit(0);
    }
    return 1;
}

void vtkMeshLess::CE(vtkPolyData* input, vtkPolyData* output) {
    output->DeepCopy(input);
    this->verPre->DeepCopy(input->GetPoints());
    vtkIdType sizePts = this->verPre->GetNumberOfPoints();
    this->verModify->DeepCopy(input->GetPoints());
    for (int i=0; i<this->ControlIds->GetNumberOfTuples(); i++) {
        double pt[3];
        this->ControlPoints->GetPoint(i, pt);
        this->verModify->SetPoint(this->ControlIds->GetValue(i), pt[0], pt[1], pt[2]);
    }
    // Get all parameters such as volume, mass, speed, positions(verPreMat and verModifyMat)
    Matrix verOriMat = pts2Matrix(this->verPre);
    Matrix verPreMat = pts2Matrix(this->verPre);
    Matrix verModifyMat = pts2Matrix(this->verModify);
    // Get the neighbor matrix, volume of bounding box, mass of each point
    Matrix neiMat = neighborMatrix(input);
    double *bounds = input->GetBounds();
    double h = (bounds[1] - bounds[0]) / this->Subdivision;
    double mass = pow((h/3),3)*5;
    // Declare matrices (will be used in the loop)
    // For the meaning of these variables, see paper "Point Based Animation" by Matthias Muller
    Matrix A(3,3);      // Sum(xij' * xij * wij)
    Matrix x(1,3);      // Sum((uj-ui) * xij * wij)
    Matrix y(1,3);      // Sum((vj-vi) * xij * wij)
    Matrix z(1,3);      // Sum((wj-wi) * xij * wij)
    Matrix xw(1,3);     // Sum(xij * wij)
    Matrix xij(1,3);    // Point i - Point j, 1-by-3 vector
    Matrix invA;        // Inverse of matrix A
    Matrix du(3,3);     // Delta U (derivative of displacement)
    Matrix J(3,3);      // Jacobian
    Matrix I(3,3);      // Identity matrix
    I.Unit();
    Matrix strain(3,3); // Strain, "epsilon" in the paper
    Matrix C(3,3);      // Approximating the constitutive law of the material
    C.Null();
    double Cfactor = this->YoungModulus / (1 - pow(this->PoissonRatio,2));
    C(0,0) = Cfactor;
    C(0,1) = Cfactor * this->PoissonRatio;
    C(1,0) = Cfactor * this->PoissonRatio;
    C(1,1) = Cfactor;
    C(2,2) = Cfactor * (1 - this->PoissonRatio) / 2;
    Matrix stress(3,3); // Stress, "delta" in the paper
    Matrix di(3,1);     // invA * (-xw)', for calculating force later
    Matrix force(3,1);
    Matrix acceleration(3,1);
    Matrix dv(3,1);     // Delta V (derivative of velocity)
    Matrix speed(sizePts, 3);
    speed.Null();
    Matrix tmpSpeed(1,3);
    Matrix tmpPos(1,3);
    // Error terms for stop criterion
    Matrix errMat = verPreMat-verModifyMat;
    double errPre = errMat.Norm() + this->Threshold;
    // Start the loop: displacement (u at t) -> derivatives (du) -> strains (epsilon)
    // -> stresses (delta) -> forces (f) -> integration (u at t+deltaT)
    //while (errPre >= errMat.Norm() && abs(errMat.normI())>this->Threshold) {
    while (1) {
        errPre = errMat.Norm();
        // Loop on every point
        for (int i=0; i<sizePts; i++) {
            A.Null();
            x.Null();
            y.Null();
            z.Null();
            xw.Null();
            double pi = 0; // Density at phyxel i
            // For every point, loop on every neighbor to get delta u (du) and wij
            for (int j=1; j<neiMat(i,0); j++) {
                // Don't consider the point itself (point i)
                if (neiMat(i,j) == i) {
                    continue;
                }
                double r = verOriMat.Distance(i,j);
                double wij = (315/64*M_PI*pow(h,9)) * pow((h*h-r*r), 3);
                pi = pi + mass * wij;
                double uj = verModifyMat(j,0) - verOriMat(j,0);
                double ui = verModifyMat(i,0) - verOriMat(i,0);
                xij(0,0) = verOriMat(j,0) - verOriMat(i,0);
                xij(0,1) = verOriMat(j,1) - verOriMat(i,1);
                xij(0,2) = verOriMat(j,2) - verOriMat(i,2);
                x = x + (uj - ui) * xij * wij;
                double vj = verModifyMat(j,1) - verOriMat(j,1);
                double vi = verModifyMat(i,1) - verOriMat(i,1);
                y = y + (vj - vi) * xij * wij;
                double wj = verModifyMat(j,2) - verOriMat(j,2);
                double wi = verModifyMat(i,2) - verOriMat(i,2);
                z = z + (wj - wi) * xij * wij;
                A = A + (~xij) * xij * wij;
                xw = xw + xij * wij;
            } // end of loop for each neighbor
            double volume = mass / pi;
            cout << pi << " " << volume << endl;
            invA = A.Inv();
            du.SetRow(0, ~(invA*(~x)));
            du.SetRow(1, ~(invA*(~y)));
            du.SetRow(2, ~(invA*(~z)));
            // Compute Jacobian, Strain, Stree
            J = du + I;
            strain = (~J)*J - I;
            stress = C * strain;
            // Force
            di = invA * (~(xw*=-1));
            force = -2 * volume * J * stress * di;
            // Integration (using euler method)
            acceleration = force / mass;
            dv = acceleration * this->DeltaT;
            tmpSpeed(0,0) = speed(i,0) + dv(0,0);
            tmpSpeed(0,1) = speed(i,1) + dv(1,0);
            tmpSpeed(0,2) = speed(i,2) + dv(2,0);
            speed.SetRow(i, tmpSpeed);
            tmpPos(0,0) = verModifyMat(i,0) + speed(i,0) * this->DeltaT;
            tmpPos(0,1) = verModifyMat(i,1) + speed(i,1) * this->DeltaT;
            tmpPos(0,2) = verModifyMat(i,2) + speed(i,2) * this->DeltaT;
            verModifyMat.SetRow(i, tmpPos);
        } // end of loop for each point
        verPreMat = verModifyMat;
        for (int i=0; i<this->ControlIds->GetNumberOfTuples(); i++) {
            double pt[3];
            this->ControlPoints->GetPoint(i, pt);
            verModifyMat(this->ControlIds->GetValue(i), 0)=pt[0];
            verModifyMat(this->ControlIds->GetValue(i), 1)=pt[1];
            verModifyMat(this->ControlIds->GetValue(i), 2)=pt[2];
        }
        errMat = verModifyMat - verPreMat;
        cout << errMat.Norm() << endl;
    } // end of loop while
    this->verModify = matrix2Pts(verModifyMat);
    output->SetPoints(this->verModify);
}

Matrix vtkMeshLess::neighborMatrix(vtkPolyData* input) {
    int ptsSize = input->GetNumberOfPoints();
    double *bounds = input->GetBounds();
    double h = (bounds[1] - bounds[0]) / this->Subdivision;
    int stepsy = ceil((bounds[3] - bounds[2]) / h);
    int stepsz = ceil((bounds[5] - bounds[4]) / h);
    int x = this->Subdivision + 2;
    int y = stepsy + 2;
    int z = stepsz + 2;
    Matrix count(x*y*z+y*z+z,1);
    count.Ones();
    for (int i=0; i<ptsSize; i++) {
        double *pt = input->GetPoint(i);
        int xt=floor((pt[0] - bounds[0]) / h) + 1;
        int yt=floor((pt[1] - bounds[2]) / h) + 1;
        int zt=floor((pt[2] - bounds[4]) / h) + 1;
        int idx = xt*(y*z) + yt*z + zt;
        count(idx,0) = count(idx,0) + 1;
    }
    int maxhash = count.normI();
    Matrix hash(x*y*z+y*z+z, maxhash);
    hash.Ones();
    for (int i=0; i<ptsSize; i++) {
        double *pt = input->GetPoint(i);
        int xt=floor((pt[0] - bounds[0]) / h) + 1;
        int yt=floor((pt[1] - bounds[2]) / h) + 1;
        int zt=floor((pt[2] - bounds[4]) / h) + 1;
        int idx = xt*(y*z) + yt*z + zt;
        int cur = hash(idx, 0);
        hash(idx, cur) = i;
        hash(idx, 0) = cur + 1;
    }
    Matrix nlistcount(ptsSize, 1);
    for (int i=0; i<hash.RowNo(); i++) {
        if (hash(i,0) == 1)
            continue;
        double ncell[19] = {i,i+1,i-1, i+1+z, i+1-z, i-1+z,
                            i-1-z, i+1+z*y, i+1-z*y, i-1+z*y, i-1-z*y,
                            i+z, i-z, i+z+z*y, i+z-z*y, i-z+z*y, i-z-z*y,
                            i+z*y, i-z*y
                           };
        for (int j=1; j<hash(i,0); j++) {
            int cur = 1;
            for (int k=0; k<19; k++) {
                if (hash(ncell[k], 0) == 1)
                    continue;
                cur = cur + hash(ncell[k], 0) - 1;
            }
            nlistcount(hash(i,j), 0) = cur;
        }
    }
    int maxnlist = nlistcount.normI();
    Matrix nlist(ptsSize, maxnlist + 1);
    nlist.Null();
    for (int i=0; i<hash.RowNo(); i++) {
        if (hash(i,0) == 1)
            continue;
        double ncell[19] = {i,i+1,i-1, i+1+z, i+1-z, i-1+z,
                            i-1-z, i+1+z*y, i+1-z*y, i-1+z*y, i-1-z*y,
                            i+z, i-z, i+z+z*y, i+z-z*y, i-z+z*y, i-z-z*y,
                            i+z*y, i-z*y
                           };
        for (int j=1; j<hash(i,0); j++) {
            int cur = 0;
            for (int k=0; k<19; k++) {
                if (hash(ncell[k], 0) == 1)
                    continue;
                for (int m=1; m<hash(ncell[k], 0); m++) {
                    nlist(hash(i,j), cur+m) = hash(ncell[k], m);
                }
                cur = cur + hash(ncell[k], 0) - 1;
            }
            nlist(hash(i,j), 0) = cur + 1;
        }
    }
    return nlist;
}

// VTK specific method:
//      This method prints the contents of this object. The method is needed
//      to conform to the VTK tools.
void vtkMeshLess::PrintSelf(ostream& os, vtkIndent indent) {
    vtkPolyDataAlgorithm::PrintSelf(os,indent);
}
