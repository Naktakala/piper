#include "ParallelPETScMatrixProxy.h"

namespace chi_math
{

ParallelPETScMatrixProxy::ParallelPETScMatrixProxy(Mat matrix,
                                         uint64_t num_local_rows,
                                         uint64_t num_local_cols,
                                         uint64_t num_global_rows,
                                         uint64_t num_global_cols,
                                         MPI_Comm comm)
  : ParallelMatrix(
      num_local_rows, num_local_cols, num_global_rows, num_global_cols, comm),
    matrix_(matrix)
{
}

void ParallelPETScMatrixProxy::Assemble(bool final)
{
  if (final)
  {
    MatAssemblyBegin(matrix_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix_, MAT_FINAL_ASSEMBLY);
  }
  else
  {
    MatAssemblyBegin(matrix_, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(matrix_, MAT_FLUSH_ASSEMBLY);
  }
}

void ParallelPETScMatrixProxy::ZeroRow(int64_t row_index, double diag_val)
{
  MatZeroRows(matrix_, 1, &row_index, diag_val, nullptr, nullptr);
}

void ParallelPETScMatrixProxy::SetValue(int64_t i, int64_t j, double value)
{
  MatSetValue(matrix_, i, j, value, INSERT_VALUES);
}

void ParallelPETScMatrixProxy::AddValue(int64_t i, int64_t j, double value)
{
  MatSetValue(matrix_, i, j, value, ADD_VALUES);
}

} // namespace chi_math