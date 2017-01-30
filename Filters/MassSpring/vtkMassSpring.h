/*********************************************************

Filename: vtkMassSpring.h

Description:

Mass Spring System. Using Runge Kutta 4 to solve this ODE:
mass*acceleration + damping*speed + spring coefficient*displacement
= external force

Input: Original object, control/fixed points id and positions

Output: deformed object

Other parameters: spring/distance coefficient, damping coefficient,
mass for all points delta t for one step in RK4, amount of steps, how
many frames needed to be saved

Usage: see python/msstest.py

Author: Shaoting Zhang

Notes: Based on modeling program of IGST in SJTU

*********************************************************/

#ifndef __vtkMassSpring_h__
#define __vtkMassSpring_h__

#include <vector>
#include "Distance.h"
#include "UGridAccelerator.h"
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkObjectFactory.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>

using namespace std;

class vtkUnstructuredGrid;

class VTK_EXPORT vtkMassSpring : public vtkUnstructuredGridAlgorithm {
public:
    static vtkMassSpring* New();
    vtkTypeMacro(vtkMassSpring,vtkUnstructuredGridAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    //get parameters for Mass Spring System
    vtkSetMacro(DistanceCoefficient, double);     // Distance constraint coefficient
    vtkSetMacro(DampingCoefficient, double);      // Damping coefficient
    vtkSetMacro(DeltaT, double);                  // dt for every step
    vtkSetMacro(Steps, int);                      // The amount of steps for simulation. Simulation time = Steps * DeltaT
    vtkSetMacro(DisplayFrames, int);              // How many key frames (for displaying later) needed. Default is 0.
    vtkSetMacro(Mass, double);                    // Mass of all points (future work: specify different masses for different "layers")

protected:
    vtkMassSpring();
    ~vtkMassSpring();

    vtkPointSet* GetControlPointsData();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int FillInputPortInformation(int, vtkInformation*);

    //data
    vtkIdTypeArray* ControlIds;
    vtkPoints* ControlPoints;
    vtkPoints* points;
    double DistanceCoefficient;
    double DeltaT;
    int Steps;
    int DisplayFrames;
    double DampingCoefficient;
    double Mass;
    double* posNext;
    double* speedNext;
    double* forces;            // array for adding up all resulting forces
    double* tk;                // for RK4
    double* tj;                // for RK4
    double* k1;                // for RK4
    double* k2;                // for RK4
    double* k3;                // for RK4
    double* k4;                // for RK4
    double* j1;                // for RK4
    double* j2;                // for RK4
    double* j3;                // for RK4
    double* j4;                // for RK4
    int numPoints;
    UGridAccelerator* gridAccel;
    vtkUnstructuredGridWriter* writer;

    // for RungeKutta4
    void Init(vtkUnstructuredGrid* input);                                     // Init data structure like posNext, speedNext, etc
    void RK4ComputeNextStep(double* posNext,     double* speedNext);           // Using Runge Kutta 4 to calculate speeds and positions of next step
    void evalSpeed(double* result, double* positions, double* speeds);         // Get speeds of all points for current step
    void evalAcceleration(double* result, double* positions, double* speeds);  // Get accelerations of all points for current step
    void addForce(vtkIdType pointId, double force[3]);                         // Add force to a specified point
    // for generate filename
    void filename(int i, char* c);    // input is the current frame "i", output is the corresponding filename "c"

private:
    vtkMassSpring(const vtkMassSpring&);            // Not implemented.
    void operator=(const vtkMassSpring&);           // Not implemented.
};

#endif

