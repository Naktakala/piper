#include "EquationSystem.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/spatial_discretization.h"

#include "chi_log.h"

#include <iomanip>

#define scint64_t static_cast<int64_t>

namespace chi_math
{
// ##################################################################
std::vector<EquationSystem::FieldBlockInfo>
EquationSystem::MakeFieldBlockInfo(const FieldList& field_list)
{
  const int num_ranks = Chi::mpi.process_count;
  const size_t num_blocks = field_list.size();

  std::vector<int64_t> blockB_n_block_local_dofs;
  std::vector<int64_t> blockB_n_block_global_dofs;

  for (const auto& field : field_list)
  {
    blockB_n_block_local_dofs.push_back(
      scint64_t(field->SDM().GetNumLocalDOFs(field->UnkManager())));
    blockB_n_block_global_dofs.push_back(
      scint64_t(field->SDM().GetNumGlobalDOFs(field->UnkManager())));
  }

  std::vector<int64_t> serial_n_block_local_dofs(num_ranks * num_blocks, 0);
  MPI_Allgather(blockB_n_block_local_dofs.data(), // sendbuf
                static_cast<int>(num_blocks),     // sendcnt
                MPI_INT64_T,                      // sendtype
                serial_n_block_local_dofs.data(), // recvbuf
                static_cast<int>(num_blocks),     // recvcount (per process)
                MPI_INT64_T,                      // recvtype
                Chi::mpi.comm);                   // communicator

  //=================================== Convert serial num_local_dofs to matrix
  typedef std::vector<int64_t> Vec64_t;
  typedef std::vector<Vec64_t> Mat64_t;

  Mat64_t mat_locP_blockB_n_block_local_dofs(num_ranks, Vec64_t(num_blocks, 0));
  size_t c = 0;
  for (int p = 0; p < num_ranks; ++p)
    for (size_t b = 0; b < num_blocks; ++b)
      mat_locP_blockB_n_block_local_dofs[p][b] = serial_n_block_local_dofs[c++];

  auto PrintMat64_t = [](const Mat64_t& mat, size_t dim_a, size_t dim_b)
  {
    std::stringstream outstr;
    for (size_t a = 0; a < dim_a; ++a)
    {
      for (size_t b = 0; b < dim_b; ++b)
      {
        outstr << "[" << a << "][" << b << "]" << std::setw(6) << mat[a][b]
               << "\n";
      }
      outstr << "\n";
    }
    return outstr.str();
  };

  if (verbosity_ >= 2)
  {
    Chi::log.Log0Warning() << "mat_locP_blockB_n_block_local_dofs:\n"
                           << PrintMat64_t(mat_locP_blockB_n_block_local_dofs,
                                           num_ranks,
                                           num_blocks)
                           << std::endl;
  }

  //=================================== Compute local offsets
  Mat64_t mat_locP_blockB_local_block_offsets(num_ranks,
                                              Vec64_t(num_blocks, 0));
  for (int p = 0; p < num_ranks; ++p)
  {
    int64_t local_offset = 0;
    for (size_t b = 0; b < num_blocks; ++b)
    {
      mat_locP_blockB_local_block_offsets[p][b] = local_offset;
      local_offset += mat_locP_blockB_n_block_local_dofs[p][b];
    }
  }

  if (verbosity_ >= 2)
  {
    Chi::log.Log0Warning() << "mat_locP_blockB_local_block_offsets:\n"
                           << PrintMat64_t(mat_locP_blockB_local_block_offsets,
                                           num_ranks,
                                           num_blocks)
                           << std::endl;
  }

  //=================================== Compute global offsets
  Mat64_t mat_locP_blockB_global_offsets(num_ranks, Vec64_t(num_blocks, 0));
  for (size_t b = 0; b < num_blocks; ++b)
  {
    int64_t block_offset = 0;
    for (int p = 0; p < num_ranks; ++p)
    {
      mat_locP_blockB_global_offsets[p][b] = block_offset;
      block_offset += mat_locP_blockB_n_block_local_dofs[p][b];
    }
  }

  if (verbosity_ >= 2)
  {
    Chi::log.Log0Warning() << "mat_locP_blockB_global_offsets:\n"
                           << PrintMat64_t(mat_locP_blockB_global_offsets,
                                           num_ranks,
                                           num_blocks)
                           << std::endl;
  }

  //=================================== Compute global offsets
  Mat64_t mat_locP_blockB_system_offsets(num_ranks, Vec64_t(num_blocks, 0));
  for (int pi = 0; pi < num_ranks; ++pi)
  {
    for (size_t bi = 0; bi < num_blocks; ++bi)
    {
      int64_t global_offset = mat_locP_blockB_local_block_offsets[pi][bi];
      for (int pj = 0; pj < pi; ++pj)
        for (size_t bj = 0; bj < num_blocks; ++bj)
          // if (bj != bi)
          global_offset += mat_locP_blockB_n_block_local_dofs[pj][bj];

      mat_locP_blockB_system_offsets[pi][bi] = global_offset;
    }
  }

  if (verbosity_ >= 2)
  {
    Chi::log.Log0Warning() << "mat_locP_blockB_system_offsets:\n"
                           << PrintMat64_t(mat_locP_blockB_system_offsets,
                                           num_ranks,
                                           num_blocks)
                           << std::endl;
  }

  //======================================= Transpose the matrix
  Mat64_t mat_locP_blockB_system_offsetsT(num_blocks, Vec64_t(num_ranks, 0));

  for (size_t b = 0; b < num_blocks; ++b)
    for (int p = 0; p < num_ranks; ++p)
    {
      mat_locP_blockB_system_offsetsT[b][p] =
        mat_locP_blockB_system_offsets[p][b];
    }

  //======================================== Determine ghost dofs
  Mat64_t blockB_ghost_ids(num_blocks);
  int64_t primary_ghost_offset = 0;
  for (size_t b = 0; b < num_blocks; ++b)
  {
    auto& field = field_list[b];
    blockB_ghost_ids[b] = field->SDM().GetGhostDOFIndices(field->UnkManager());
    primary_ghost_offset += scint64_t(blockB_n_block_local_dofs[b]);
  }
  int64_t secondary_ghost_offset = 0;
  Vec64_t blockB_ghost_local_offsets(num_blocks, 0);
  for (size_t b = 0; b < num_blocks; ++b)
  {
    blockB_ghost_local_offsets[b] =
      secondary_ghost_offset + primary_ghost_offset;
    secondary_ghost_offset += scint64_t(blockB_ghost_ids[b].size());
  }

  std::vector<EquationSystem::FieldBlockInfo> field_block_info_list;

  const int location_id = Chi::mpi.location_id;
  for (size_t b = 0; b < num_blocks; ++b)
  {
    auto& field = field_list[b];
    const int64_t num_local_dofs = blockB_n_block_local_dofs[b];
    const int64_t num_global_dofs = blockB_n_block_global_dofs[b];

    std::vector<std::pair<int64_t, int64_t>> locP_block_span(num_ranks, {0, 0});
    for (int p = 0; p < num_ranks; ++p)
    {
      const int64_t begin = mat_locP_blockB_global_offsets[p][b];
      const int64_t end = begin + mat_locP_blockB_n_block_local_dofs[p][b];
      locP_block_span[p] = {begin, end};
    }

    field_block_info_list.push_back(std::move(
      FieldBlockInfo{field,
                     num_local_dofs,
                     num_global_dofs,
                     mat_locP_blockB_local_block_offsets[location_id][b],
                     blockB_ghost_local_offsets[b],
                     blockB_ghost_ids[b],
                     locP_block_span,
                     mat_locP_blockB_system_offsetsT[b]}));
  }

  return field_block_info_list;
}

int64_t EquationSystem::MapBlockGlobalIDToSystem(const FieldBlockInfo& info,
                                                 int64_t block_global_id) const
{
  const size_t num_ranks = info.locP_global_block_span_.size();

  const auto& span = info.locP_global_block_span_;
  const auto& system_offset = info.locP_system_offsets_;

  int64_t dof_map = -1;
  for (size_t p = 0; p < num_ranks; ++p)
    if (block_global_id >= span[p].first and block_global_id < span[p].second)
    {
      dof_map = block_global_id - span[p].first + system_offset[p];
      break;
    }

  ChiLogicalErrorIf(dof_map >= num_globl_dofs_,
                    "dof_map " + std::to_string(dof_map) + " >= global-size " +
                      std::to_string(num_globl_dofs_));

  if (dof_map >= 0) return dof_map;

  for (size_t p = 0; p < num_ranks; ++p)
  {
    Chi::log.LogAllError() << "p=" << p << " begin=" << span[p].first
                           << " end=" << span[p].second
                           << " block_gid=" << block_global_id;
  }

  ChiLogicalError("Failed to map block global-id to system global-id");
}

int64_t EquationSystem::MapBlockLocalIDToSystem(const FieldBlockInfo& info,
                                                int64_t block_local_id)
{
  int64_t dof_map = 0;
  if (block_local_id < info.num_local_dofs_)
    dof_map = block_local_id + info.local_offset_;
  else // its ghost local id
  {
    const int64_t ghost_sub_id = block_local_id - info.num_local_dofs_;
    dof_map = ghost_sub_id + info.ghost_local_offset_;
  }

  ChiLogicalErrorIf(dof_map < 0,
                    "dof_map < zero for block_local_id " +
                      std::to_string(block_local_id));
  const int64_t num_local_values =
    scint64_t(main_solution_vector_->RawValues().size());
  ChiLogicalErrorIf(dof_map >= num_local_values,
                    "dof_map " + std::to_string(dof_map) +
                      " >= num_local_values " +
                      std::to_string(num_local_values));

  return dof_map;
}

} // namespace chi_math