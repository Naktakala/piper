#include "ParallelMatrix.h"

namespace chi_math
{

ParallelMatrix::ParallelMatrix(uint64_t num_local_rows,
                               uint64_t num_local_cols,
                               uint64_t num_global_rows,
                               uint64_t num_global_cols,
                               MPI_Comm comm)
  : num_local_rows_(num_local_rows),
    num_local_cols_(num_local_cols),
    num_global_rows_(num_global_rows),
    num_global_cols_(num_global_cols),
    comm_(comm)
{
}

uint64_t ParallelMatrix::NumLocalRows() const { return num_local_rows_; }
uint64_t ParallelMatrix::NumGlobalRows() const { return num_global_rows_; }
uint64_t ParallelMatrix::NumLocalColumns() const { return num_local_cols_; }
uint64_t ParallelMatrix::NumGlobalColumns() const { return num_global_cols_; }

MPI_Comm ParallelMatrix::Comm() const { return comm_; }

void ParallelMatrix::SetSparsityPattern(
  const std::function<void()>& sparsity_function)
{
  sparsity_function();
}

} // namespace chi_math