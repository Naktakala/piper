#ifndef CHITECH_PARALLELPETSCMATRIXPROXY_H
#define CHITECH_PARALLELPETSCMATRIXPROXY_H

#include "ParallelMatrix.h"

#include <petsc.h>

namespace chi_math
{

class ParallelPETScMatrixProxy : public ParallelMatrix
{
public:
  ParallelPETScMatrixProxy(Mat matrix,
                      uint64_t num_local_rows,
                      uint64_t num_local_cols,
                      uint64_t num_global_rows,
                      uint64_t num_global_cols,
                      MPI_Comm comm);

  void Assemble() override;

  void SetValue(int64_t i, int64_t j, double value) override;
  void AddValue(int64_t i, int64_t j, double value) override;

protected:
  Mat matrix_;
};

}

#endif // CHITECH_PARALLELPETSCMATRIXPROXY_H
