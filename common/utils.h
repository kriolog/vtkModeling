#ifndef __utils_h__
#define __utils_h__

#include <vtkIdTypeArray.h>

vtkIdType IsId(vtkIdTypeArray* ids, vtkIdType id)
{
  for(vtkIdType i = 0, N = ids->GetNumberOfTuples(); i < N; ++i) {
    if (id == ids->GetValue(i))
      return i;
  }
  return (-1);
}

#endif
