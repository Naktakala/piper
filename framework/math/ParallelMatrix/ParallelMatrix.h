#ifndef CHITECH_PARALLELMATRIX_H
#define CHITECH_PARALLELMATRIX_H

#include "chi_mpi.h"

#include <cstdint>
#include <functional>

namespace chi_math
{

struct ParallelMatrixSparsityPattern
{
  std::vector<int64_t> nodal_nnz_in_diag_;
  std::vector<int64_t> nodal_nnz_off_diag_;
};

class ParallelMatrix
{
public:
  uint64_t NumLocalRows() const;
  uint64_t NumGlobalRows() const;
  uint64_t NumLocalColumns() const;
  uint64_t NumGlobalColumns() const;

  MPI_Comm Comm() const;

  virtual void Assemble(bool final) = 0;
  virtual void ZeroRow(int64_t row_index, double diag_val) = 0;

  virtual void SetValue(int64_t i, int64_t j, double value) = 0;
  virtual void AddValue(int64_t i, int64_t j, double value) = 0;

  virtual void
  SetSparsityPattern(const std::function<void()>& sparsity_function);

  virtual ~ParallelMatrix() = default;

protected:
  ParallelMatrix(uint64_t num_local_rows,
                 uint64_t num_local_cols,
                 uint64_t num_global_rows,
                 uint64_t num_global_cols,
                 MPI_Comm comm);

  const uint64_t num_local_rows_;
  const uint64_t num_local_cols_;
  const uint64_t num_global_rows_;
  const uint64_t num_global_cols_;
  const MPI_Comm comm_;
};

} // namespace chi_math

#endif // CHITECH_PARALLELMATRIX_H
