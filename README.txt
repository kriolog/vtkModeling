INSTALL

Step 0:
Basically you need to install cmake, MSVC or GCC + Make, and configure
VTK5.0. If you wanna run the test easily, you can also install Python.

Step 1:
Linux:
cd vtkExtend/build
ccmake ../
# change the "VTK_DIR" to your vtk build directory
# Turn off "WRAP_PYTHON" if you don't have Python
# Press c
# Press g
# quit CMake mode
make

MSVC:
CMake the project, set the build dir to vtkExtend/build
Turn off "BUILD_SHARED_LIBS" (otherwise may cause errors)
Turn off "WRAP_PYTHON" if you wanna compile with debug mode (otherwise
Python wrap tools may report errors)
Generate project
compile with release mode if enable "WRAP_PYTHON"

Step2: demo
Python demo:
python/*test.py

C++ demo:
*/test_*.cpp

-------------------------------------------------------------

Author: Shaoting Zhang
Contact method: rutgers.shaoting@gmail.com
Homepage in Rutgers: paul.rutgers.edu/~shaoting


-------------------------------------------------------------

Cite this work:

S. Zhang, X. Wang, D. Metaxas, T. Chen and L. Axel, 
LV Surface Reconstruction From Sparse tMRI Using Laplacian Surface 
Deformation And Optimization, In Proc. of ISBI'09, Boston, US, June 28 - July 1, 2009.

-------------------------------------------------------------

How to extend VTK (in case you're interested):

		VTKModeling, an example of extending VTK

During the summer I implemented a modeling toolkit by extending
Visualization ToolKit. Basically each algorithm works like a filter.
It gets a dataset and bunches of parameters as input, and generate a
deformed dataset as the output. It's very easy to use, so I also made
a few applications. All source codes can be downloaded in
Sourceforge.net. They're under GPL, so feel free to try and add new
stuff. Here I wrote down the whole extending procedure.

0.
Install CMake, VTK5.0, Python (if you need), MSVC or Make + GCC + GDB

1.
Figure out the proper input and output for your algorithm. In my case,
vtkPolyData is preferable for laplacian surface editing since I need
the connectivity information. vtkUnstructuredGrid is better for
volume-reserving mass spring system. vtkPoints is fine for meshless
model.

2.
Extend proper subclass of vtkAlgorithm according to the input and
output, such as vtkPolyDataAlgorithm for laplacian surface editing,
vtkUnstructuredGridAlgorithm for mass spring system. Of course you can
always change the data type of input and output by some extra codes,
but sticking to correct class can save time.

3.
Write your own codes based on the skeleton one (or demo codes, like
all subclasses of vtkAlgorithm). Check out this example firstly
(vtkLaplacianSurface):
(vtkLaplacianSurface.h)
#include <vtkPolyDataAlgorithm.h>
class vtkPolyData;
class VTK_EXPORT vtkLaplacianSurface : public vtkPolyDataAlgorithm {
public:
    static vtkLaplacianSurface* New();
    vtkTypeMacro(vtkLaplacianSurface,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    //Get parameters for lse. This declaration is enough, no definition needed.
    vtkSetObjectMacro(ControlIds, vtkIdList);    // Input, control points' id
    vtkSetObjectMacro(ControlPoints, vtkPoints);// Input, control points' coordinates
protected:
    vtkLaplacianSurface();
    ~vtkLaplacianSurface();
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int FillInputPortInformation(int, vtkInformation*);
    //data
    vtkIdList* ControlIds;   // control points' id
    vtkPoints* ControlPoints;// control points' coordinates
private:
    vtkLaplacianSurface(const vtkLaplacianSurface&);    // Not implemented.
    void operator=(const vtkLaplacianSurface&);         // Not implemented.
};
(vtkLaplacianSurface.cpp), need to define one important method:
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
    doLSE(input, output);
    return 1;
}

Here "input" holds the pointer to the input dataset. We CANNOT change
it, but we can always make a "deepcopy" and deal with that copy one.
"output" holds the pointer to the output dataset of this filter, which
may be used by other filter later. It's empty in the beginning. We can
feed in dataset at anytime.

4.
Write CMakelist to maintain project conveniently. See this example
(for vtkLaplacianSurface):

(CMakeList.txt)
SET (LIBRARY_NAME
    lse_vtk
)
SET (SOURCES
    vtkLaplacianSurface.cpp
    MathTools.cpp
)
SET_SOURCE_FILES_PROPERTIES(
    MathTools.cpp
  PROPERTIES
    WRAP_EXCLUDE 1
)
SET (WRAPPED_LIBRARIES
    vtkCommon
    vtkFiltering
    vtkVolumeRendering
)
SET (LIBRARIES
    ${WRAPPED_LIBRARIES}
)
SET (INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/${LIBRARY_NAME}
    ${PROJECT_SOURCE_DIR}/include
)
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
ADD_EXECUTABLE(test_LS test_LS.cpp)
TARGET_LINK_LIBRARIES(test_LS vtkRendering lse_vtk)
INCLUDE(../CMakeCommon.txt)

Two things worth mention: "WRAP_EXCLUDE 1" means that "don't wrap
these files for script languages like python". Normally these files
are not satisfy warping criterion, besides, they just like private
methods and useless for python. Of course they can be invoked by other
C++ files, like vtkLaplacianSurface.cpp in this case. "ADD_EXECUTABLE"
will generate a new C++ project. It's helpful for debugging.

5. 
CMake your codes, which will generate a makefile in Linux or a MSVC
project in Windows. Compile and run it. Everything is ok in Linux, but
I got one weird error in MSVC: Python Wrap Tool will return fatal
errors when compiling in Debug mode. Release mode is ok. So, to debug
in MSVC, we got to disable Python wrap. It's easy to do. Just add
several lines to CMakeList.txt: OPTION(WRAP_PYTHON "WRAP_PYTHON."
${WRAP_PYTHON})

That's it.

---------------------------- More details ----------------------------

Project link:
http://sourceforge.net/projects/vtkextend/

Deformation algorithms implemented so far:
Laplacian Surface Editing and Optimization
Meshless model with LSE, or with Continuum Elasticity
Moving Least Square
Mass Spring System (almost the same as the Robust Deformable Model in
IGST/SJTU)

Matrix calculation:
Currently I'm using "Matrix TCL Lite". It's easy to use (only one .h
file), and works like Matlab in some way. But it's really slow for big
matrix operations. We applied these algorithms on clinical data with
around 3,000 vertices. It takes hours for laplacian surface
optimization and 10 minutes for editing. Matrix TCL is for dense
matrix only, without any optimization for sparse one. That's why it's
also slow for our sparse matrix calculation. Some libraries are
available for sparse matrix. They're awesome, extremely high
efficiency, but really hard to configure and use. Finally I decided to
write several functions for sparse matrix calculation. Basically I'm
using conjugate gradient for solving Ax=b, where A is a sparse
symmetric definite matrix, and using a simple idea to solve A*B, where
A is a sparse matrix. So we can always change Ax=b to A^T*A*x = A^T*b,
where A is sparse (and A^T*A will be sparse symmetric definite). I put
them in Matrix TCL, named as SparseMul and SSDSolver. The running
decreases from hours to one minute.

Debug with GDB:
Three problems so far, only one solved.
1. Operator overload: "p aMatrix(1,1)" won't work if "(" and ")" are
overloaded operator. Try "p aMatrix.'Operator()'(1,1)" instead.
2. Inline function: due to GCC's optimization, inline functions won't
be linked. We cannot invoke them in GDB. I'm using CMake, no idea how
to control the optimization.
3. Methods in VTK: cannot invoke VTK methods in GDB through object
reference. How come... Thanks a lot for suggestions.

Implementation limitations:
1. Moving Least Square only works for 2D (cannot figure out the
transformation for 3D case).
2. Meshless with Continuum Elasticity may not be correct.
3. Added several functions in Matrix TCL, but didn't write exception handling.

References:
Olga Sorkine, "Laplacian Surface Editing"
Andrew Nealen, "Laplacian Surface Optimization"
Xiaoxu Wang, "LV Motion and Strain Computation from tMRI based on
Meshless Deformable Models"
Matthias Muller, "Point Based Animation of Elastic, Plastic and Melting Object"
Scott Schaefer, "Image Deformation Using Moving Least Square"
Jonathan Richard Shewchuk, "An Introduction to the Conjugate Gradient
Method Without the Agonizing Pain"
Mark Meyery and Haeyoung Leez, "Generalized Barycentric Coordinates on
Irregular Polygons"

