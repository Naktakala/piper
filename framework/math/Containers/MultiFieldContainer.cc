#include "MultiFieldContainer.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/spatial_discretization.h"
#include "math/ParallelMatrix/ParallelMatrix.h"

#include "ChiObjectFactory.h"

#include "chi_log.h"

#define scint64_t static_cast<int64_t>
#define cint64_t const int64_t

namespace chi_math
{

RegisterChiObject(chi_math, MultiFieldContainer);

chi::InputParameters MultiFieldContainer::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription(
    "A container to hold one or more grid-based field functions. The fields"
    "can have differing spatial discretizations. The container also allows for"
    "easy definition of vectors and block vectors.");

  params.AddOptionalParameter("verbosity", 0, "Verbosity level.");

  params.AddRequiredParameterArray(
    "fields", "An array of FieldFunctionGridBased handles or names.");

  return params;
}

MultiFieldContainer::MultiFieldContainer(const chi::InputParameters& params)
  : ChiObject(params),
    verbosity_(params.GetParamValue<int>("verbosity")),
    field_block_info_(std::move(MakeFieldBlockInfo(
      BuildFieldList(params.GetParam("fields")), verbosity_))),
    total_num_local_dofs_(DetermineNumLocalDofs(field_block_info_)),
    total_num_global_dofs_(DetermineNumGlobalDofs(field_block_info_)),
    system_ghost_ids_(BuildSystemGhostIDs())
{
}

// ##################################################################
/**Makes a list of grid-base field functions given either an array of
 * warehouse handles or field-function names.*/
FieldList
MultiFieldContainer::BuildFieldList(const chi::ParameterBlock& param_array)
{
  ChiInvalidArgumentIf(
    param_array.Type() != chi::ParameterBlockType::ARRAY,
    "\"fields\" parameter for EquationSystem must be of type ARRAY.");

  ChiInvalidArgumentIf(
    param_array.NumParameters() == 0,
    "\"fields\" parameter for EquationSystem must not be empty.");

  std::vector<std::string> field_names_processed;
  auto CheckDuplicateName = [&field_names_processed](const std::string& ff_name)
  {
    return std::any_of(field_names_processed.begin(),
                       field_names_processed.end(),
                       [&ff_name](const std::string& name)
                       { return name == ff_name; });
  };

  const auto values_type = param_array.GetParam(0).Type();
  if (values_type == chi::ParameterBlockType::INTEGER)
  {
    FieldList field_list;
    for (const auto& index_param : param_array)
    {
      const size_t handle = index_param.GetValue<size_t>();
      auto field_base_ptr =
        Chi::GetStackItemPtr(Chi::field_function_stack, handle, __FUNCTION__);

      auto field_ptr =
        std::dynamic_pointer_cast<chi_physics::FieldFunctionGridBased>(
          field_base_ptr);

      ChiLogicalErrorIf(not field_ptr,
                        "Field function \"" + field_base_ptr->TextName() +
                          "\" is not of required type FieldFunctionGridBased.");

      ChiInvalidArgumentIf(CheckDuplicateName(field_ptr->TextName()),
                           "Duplicate field \"" + field_ptr->TextName() +
                             "\" encountered.");

      field_list.push_back(field_ptr);
      field_names_processed.push_back(field_ptr->TextName());
    }

    return field_list;
  }
  else if (values_type == chi::ParameterBlockType::STRING)
  {
    FieldList field_list;
    for (const auto& index_param : param_array)
    {
      const std::string ff_name = index_param.GetValue<std::string>();

      std::shared_ptr<chi_physics::FieldFunctionGridBased> field_ptr = nullptr;

      for (auto& ff : Chi::field_function_stack)
        if (ff->TextName() == ff_name)
        {
          field_ptr =
            std::dynamic_pointer_cast<chi_physics::FieldFunctionGridBased>(ff);
          ChiLogicalErrorIf(
            not field_ptr,
            "Field function \"" + ff_name +
              "\" is not of required type FieldFunctionGridBased.");
        }

      ChiInvalidArgumentIf(not field_ptr,
                           "Field function \"" + ff_name +
                             "\" was not found in field-function warehouse.");

      ChiInvalidArgumentIf(CheckDuplicateName(field_ptr->TextName()),
                           "Duplicate field \"" + field_ptr->TextName() +
                             "\" encountered.");

      field_list.push_back(field_ptr);
      field_names_processed.push_back(field_ptr->TextName());
    }

    return field_list;
  }
  else
    ChiInvalidArgument(
      "The elements of the \"variables\" parameter for "
      "EquationSystem can only be INTEGER or STRING. This is either a handle "
      "to a FieldFunctionGridBased or a name to an existing one.");
}

// ##################################################################
std::vector<MultiFieldContainer::FieldBlockInfo>
MultiFieldContainer::MakeFieldBlockInfo(const FieldList& field_list,
                                        int verbosity)
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

  if (verbosity >= 2)
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

  if (verbosity >= 2)
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

  if (verbosity >= 2)
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

  if (verbosity >= 2)
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

  std::vector<FieldBlockInfo> field_block_info_list;

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

// ##################################################################
int64_t MultiFieldContainer::DetermineNumLocalDofs(
  const std::vector<FieldBlockInfo>& field_block_info)
{
  size_t num_dofs = 0;
  for (const auto& field_block : field_block_info)
    num_dofs += field_block.num_local_dofs_;

  return static_cast<int64_t>(num_dofs);
}
// ##################################################################
int64_t MultiFieldContainer::DetermineNumGlobalDofs(
  const std::vector<FieldBlockInfo>& field_block_info)
{
  size_t num_dofs = 0;
  for (const auto& field_block : field_block_info)
    num_dofs += field_block.num_global_dofs_;

  return static_cast<int64_t>(num_dofs);
}

// ##################################################################
const std::vector<int64_t>& MultiFieldContainer::GetSystemGhostIDs() const
{
  return system_ghost_ids_;
}

// ##################################################################
/**Uses the underlying system to build a sparsity pattern.*/
ParallelMatrixSparsityPattern MultiFieldContainer::BuildMatrixSparsityPattern() const
{
  std::vector<int64_t> master_row_nnz_in_diag;
  std::vector<int64_t> mastero_row_nnz_off_diag;

  auto& a1 = master_row_nnz_in_diag;
  auto& a2 = mastero_row_nnz_off_diag;
  for (const auto& field_info : field_block_info_)
  {
    auto& field = field_info.field_;
    auto& sdm = field->SDM();
    auto& uk_man = field->UnkManager();

    std::vector<int64_t> nodal_nnz_in_diag;
    std::vector<int64_t> nodal_nnz_off_diag;
    sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, uk_man);

    // Now append these to the master list
    auto& b1 = nodal_nnz_in_diag;
    auto& b2 = nodal_nnz_off_diag;

    using std::begin, std::end;
    a1.insert(end(a1), begin(b1), end(b1));
    a2.insert(end(a2), begin(b2), end(b2));
  }

  return {master_row_nnz_in_diag, mastero_row_nnz_off_diag};
}

// ##################################################################
const MultiFieldContainer::FieldBlockInfo&
MultiFieldContainer::GetFieldBlockInfo(size_t field_id) const
{
  return field_block_info_.at(field_id);
}

// ##################################################################
const chi_mesh::MeshContinuum& MultiFieldContainer::GetSystemCommonGrid() const
{
  return field_block_info_.front().field_->SDM().Grid();
}

// ##################################################################
int64_t
MultiFieldContainer::MapFieldGlobalIDToSystem(const FieldBlockInfo& info,
                                              int64_t field_global_id) const
{
  const size_t num_ranks = info.locP_global_block_span_.size();

  const auto& span = info.locP_global_block_span_;
  const auto& system_offset = info.locP_system_offsets_;

  int64_t dof_map = -1;
  for (size_t p = 0; p < num_ranks; ++p)
    if (field_global_id >= span[p].first and field_global_id < span[p].second)
    {
      dof_map = field_global_id - span[p].first + system_offset[p];
      break;
    }

  ChiLogicalErrorIf(dof_map >= total_num_global_dofs_,
                    "dof_map " + std::to_string(dof_map) + " >= global-size " +
                      std::to_string(total_num_global_dofs_));

  if (dof_map >= 0) return dof_map;

  for (size_t p = 0; p < num_ranks; ++p)
  {
    Chi::log.LogAllError() << "p=" << p << " begin=" << span[p].first
                           << " end=" << span[p].second
                           << " field_gid=" << field_global_id;
  }

  ChiLogicalError("Failed to map field global-id to system global-id");
}

// ##################################################################
int64_t
MultiFieldContainer::MapFieldLocalIDToSystem(const FieldBlockInfo& info,
                                             int64_t field_local_id) const
{
  int64_t dof_map;
  if (field_local_id < info.num_local_dofs_)
    dof_map = field_local_id + info.local_offset_;
  else // its ghost local id
  {
    const int64_t ghost_sub_id = field_local_id - info.num_local_dofs_;
    dof_map = ghost_sub_id + info.ghost_local_offset_;
  }

  ChiLogicalErrorIf(dof_map < 0,
                    "dof_map < zero for block_local_id " +
                      std::to_string(field_local_id));
  const int64_t num_local_values =
    scint64_t(total_num_local_dofs_ + system_ghost_ids_.size());
  ChiLogicalErrorIf(dof_map >= num_local_values,
                    "dof_map " + std::to_string(dof_map) +
                      " >= num_local_values " +
                      std::to_string(num_local_values));

  return dof_map;
}

// ##################################################################
std::vector<int64_t> MultiFieldContainer::BuildSystemGhostIDs()
{
  std::vector<int64_t> ghost_ids;

  for (const auto& info : field_block_info_)
    for (const int64_t field_ghost_id : info.ghost_ids_)
      ghost_ids.push_back(MapFieldGlobalIDToSystem(info, field_ghost_id));

  return ghost_ids;
}

// ##################################################################
int64_t MultiFieldContainer::MapDOF(const chi_mesh::Cell& cell,
                                    size_t i,
                                    const FieldBlockInfo& info,
                                    size_t component_id) const
{
  cint64_t field_dof_id = info.field_->SDM().MapDOF(
    cell, i, info.field_->UnkManager(), 0, component_id);

  cint64_t dof_id = MapFieldGlobalIDToSystem(info, field_dof_id);

  return dof_id;
}

// ##################################################################
int64_t MultiFieldContainer::MapDOFLocal(const chi_mesh::Cell& cell,
                                         size_t i,
                                         const FieldBlockInfo& info,
                                         size_t component_id) const
{
  cint64_t field_dof_id = info.field_->SDM().MapDOFLocal(
    cell, i, info.field_->UnkManager(), 0, component_id);

  cint64_t dof_id = MapFieldLocalIDToSystem(info, field_dof_id);

  return dof_id;
}

} // namespace chi_math