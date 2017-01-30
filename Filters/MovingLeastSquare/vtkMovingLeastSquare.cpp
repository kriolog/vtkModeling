#include "vtkMovingLeastSquare.h"

#include <math.h>
#include <iostream>
#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPoints.h>
#include <utils.h>

using namespace std;

vtkStandardNewMacro(vtkMovingLeastSquare);


// Default constructor.
vtkMovingLeastSquare::vtkMovingLeastSquare() {
    this->P = vtkPoints::New();
    this->pts = vtkPoints::New();
    this->IsAffine = 0;
    this->IsSimilarity = 0;
    this->IsRigid = 0;
    this->Dimension = 2;
    this->Alpha = 1.5;

    this->SetNumberOfInputPorts(2);
}

// Destrucor
vtkMovingLeastSquare::~vtkMovingLeastSquare() {
    this->P->Delete();
    this->pts->Delete();
    delete PMat;
    delete QMat;
    delete ptMat;
}

vtkPointSet* vtkMovingLeastSquare::GetControlPointsData()
{
  if (this->GetNumberOfInputConnections(1) < 1) {
    return NULL;
  }

  return vtkPointSet::SafeDownCast(this->GetInputDataObject(1, 0));
}

// VTK specific method:
//      This method is called when the pipeline is calculated.
int vtkMovingLeastSquare::RequestData(
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
    this->Q = ControlPointSet->GetPoints();
    this->Ids
      = vtkIdTypeArray::SafeDownCast(this->GetInputArrayToProcess(0, ControlPointSet));

    output->DeepCopy(input);
    this->pts = output->GetPoints();
    for (int i=0; i<this->Ids->GetNumberOfTuples(); i++) {
        this->P->InsertNextPoint(pts->GetPoint(this->Ids->GetValue(i)));
    }
    PMat = pts2Matrix(this->P);
    QMat = pts2Matrix(this->Q);
    ptMat = pts2Matrix(this->pts);

    for (int i=0; i<this->pts->GetNumberOfPoints(); i++) {
        if (IsId(this->Ids, i) != -1) {
            (*ptMat)(i,0)=(*QMat)(IsId(this->Ids, i), 0);
            (*ptMat)(i,1)=(*QMat)(IsId(this->Ids, i), 1);
            if (this->Dimension == 3) {
                (*ptMat)(i,2)=(*QMat)(IsId(this->Ids, i), 2);
            }
            continue;
        }
        Matrix ptTemp(1, this->Dimension);
        ptTemp.Null();
        if (this->IsAffine == 1) {
            ptTemp = affine(ptMat->GetRow(i), PMat, QMat);
        } else if (this->IsSimilarity == 1) {
            ptTemp = similarity(ptMat->GetRow(i), PMat, QMat);
        } else if (this->IsRigid == 1) {
            ptTemp = rigid(ptMat->GetRow(i), PMat, QMat);
        }
        (*ptMat)(i,0)=ptTemp(0, 0);
        (*ptMat)(i,1)=ptTemp(0, 1);
        if (this->Dimension == 3) {
            (*ptMat)(i,2)=ptTemp(0, 2);
        }
    }
    this->pts = matrix2Pts(ptMat);
    output->SetPoints(this->pts);
    return 1;
}

Matrix vtkMovingLeastSquare::affine(Matrix v, Matrix *p, Matrix *q) {
    Matrix pStar(1,this->Dimension);
    Matrix qStar(1,this->Dimension);
    Matrix fa(1,this->Dimension);
    pStar.Null();
    qStar.Null();
    fa.Null();
    double wSum = 0;
    for (unsigned int i=0; i<p->RowNo(); i++) {
        double w = weight(p->GetRow(i), v);
        pStar = pStar + p->GetRow(i) * w;
        qStar = qStar + q->GetRow(i) * w;
        wSum = wSum + w;
    }
    pStar = pStar/wSum;
    qStar = qStar/wSum;

    for (unsigned int j=0; j<p->RowNo(); j++) {
        Matrix pwp(this->Dimension,this->Dimension);
        pwp.Null();
        for (unsigned int i=0; i<p->RowNo(); i++) {
            Matrix piHat = p->GetRow(i) - pStar;
            double wi = weight(p->GetRow(i), v);
            pwp = pwp + wi*(~piHat)*piHat;
        }
        Matrix pjHat = p->GetRow(j) - pStar;
        Matrix Aj = (v-pStar)*(pwp.Inv())*(~pjHat);
        Matrix qjHat = q->GetRow(j) - qStar;
        fa = fa + Aj*qjHat;
    }
    fa = fa+qStar;
    return fa;
}

Matrix vtkMovingLeastSquare::similarity(Matrix v, Matrix *p, Matrix *q) {
    Matrix pStar(1,this->Dimension);
    Matrix qStar(1,this->Dimension);
    Matrix fs(1,this->Dimension);
    pStar.Null();
    qStar.Null();
    fs.Null();
    double wSum = 0;
    double us = 0;
    for (unsigned int i=0; i<p->RowNo(); i++) {
        double w = weight(p->GetRow(i), v);
        pStar = pStar + p->GetRow(i) * w;
        qStar = qStar + q->GetRow(i) * w;
        wSum = wSum + w;
    }
    pStar = pStar/wSum;
    qStar = qStar/wSum;

    for (unsigned int i=0; i<p->RowNo(); i++) {
        double wi = weight(p->GetRow(i), v);
        Matrix piHat = p->GetRow(i) - pStar;
        Matrix ppt = piHat*(~piHat);
        us = us + wi*ppt(0,0);
    }

    for (unsigned int i=0; i<p->RowNo(); i++) {
        double wi = weight(p->GetRow(i), v);
        Matrix piHat = p->GetRow(i) - pStar;
        Matrix qiHat = q->GetRow(i) - qStar;
        Matrix vp = v - qStar;
        Matrix Ai = wi*simOper(piHat)*(~(simOper(vp)));
        fs = fs + qiHat * (1/us * Ai);
    }
    fs = fs + qStar;
    return fs;
}

Matrix vtkMovingLeastSquare::rigid(Matrix v, Matrix *p, Matrix *q) {
    Matrix pStar(1,this->Dimension);
    Matrix qStar(1,this->Dimension);
    Matrix fr(1,this->Dimension);
    Matrix frVec(1, this->Dimension);
    pStar.Null();
    qStar.Null();
    fr.Null();
    frVec.Null();
    double wSum = 0;
    for (unsigned int i=0; i<p->RowNo(); i++) {
        double w = weight(p->GetRow(i), v);
        pStar = pStar + p->GetRow(i) * w;
        qStar = qStar + q->GetRow(i) * w;
        wSum = wSum + w;
    }
    pStar = pStar/wSum;
    qStar = qStar/wSum;

    for (unsigned int i=1; i<p->RowNo(); i++) {
        double wi = weight(p->GetRow(i), v);
        Matrix piHat = p->GetRow(i) - pStar;
        Matrix qiHat = q->GetRow(i) - qStar;
        Matrix vp = v - pStar;
        Matrix Ai = wi*simOper(piHat)*(~(simOper(vp)));
        frVec = frVec + qiHat*Ai;
    }
    Matrix vp = v - pStar;
    fr = vp.Norm() * frVec / frVec.Norm() + qStar;
    return fr;
}

double vtkMovingLeastSquare::weight(Matrix pi, Matrix v) {
    Matrix piv = pi - v;
    double w = 1 / pow(piv.Norm(), 2*this->Alpha);
    return w;
}

Matrix vtkMovingLeastSquare::simOper(Matrix m) {
    Matrix s(this->Dimension, this->Dimension);
    s.Null();
    if (this->Dimension == 2) {
        s(0,0) = m(0,0);
        s(0,1) = m(0,1);
        s(1,0) = m(0,1);
        s(1,1) = -m(0,0);
    }
    return s;
}

vtkPoints* vtkMovingLeastSquare::matrix2Pts(Matrix* m) {
    vtkPoints* pts = vtkPoints::New();
    if (this->Dimension == 3) {
        for (unsigned i=0; i<m->RowNo(); i++) {
            pts->InsertNextPoint((*m)(i,0), (*m)(i,1), (*m)(i,2));
        }
    } else if (this->Dimension == 2) {
        for (unsigned i=0; i<m->RowNo(); i++) {
            pts->InsertNextPoint((*m)(i,0), (*m)(i,1), 0);
        }
    }
    return pts;
}

Matrix* vtkMovingLeastSquare::pts2Matrix(vtkPoints* v) {
    vtkIdType sizePts = v->GetNumberOfPoints();
    Matrix* vMat = new Matrix(sizePts,this->Dimension);
    vMat->Null();
    for (int i=0; i<sizePts; i++) {
        double pt[3];
        v->GetPoint(i, pt);
        if (this->Dimension == 3) {
            (*vMat)(i,0)=pt[0];
            (*vMat)(i,1)=pt[1];
            (*vMat)(i,2)=pt[2];
        } else if (this->Dimension == 2) {
            (*vMat)(i,0)=pt[0];
            (*vMat)(i,1)=pt[1];
        }
    }
    return vMat;
}

// VTK specific method:
//      This method sets the input ports for this pipeline element.
int vtkMovingLeastSquare::FillInputPortInformation(int port, vtkInformation* info) {
  if (port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  else
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");

  return 1;
}


// VTK specific method:
//      This method prints the contents of this object. The method is needed
//      to conform to the VTK tools.
void vtkMovingLeastSquare::PrintSelf(ostream& os, vtkIndent indent) {
    vtkPolyDataAlgorithm::PrintSelf(os,indent);
}
