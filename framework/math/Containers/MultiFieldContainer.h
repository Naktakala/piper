#ifndef CHITECH_MULTIFIELDCONTAINER_H
#define CHITECH_MULTIFIELDCONTAINER_H

#include "ChiObject.h"

namespace chi_physics
{
class FieldFunctionGridBased;
}

namespace chi_mesh
{
class Cell;
class MeshContinuum;
} // namespace chi_mesh

namespace chi_math
{

typedef std::vector<std::shared_ptr<chi_physics::FieldFunctionGridBased>>
  FieldList;

/**A container to hold one or more grid-based field functions. The fields
can have differing spatial discretizations. The container also allows for
easy definition of vectors and block vectors.*/
class MultifieldContainer : public ChiObject
{
public:
  struct FieldBlockInfo;
  static chi::InputParameters GetInputParameters();
  explicit MultifieldContainer(const chi::InputParameters& params);

  int64_t TotalNumLocalDOFs() const { return total_num_local_dofs_; }
  int64_t TotalNumGlobalDOFs() const { return total_num_global_dofs_; }

  const std::vector<int64_t>& GetSystemGhostIDs() const;

  const FieldBlockInfo& GetFieldBlockInfo(size_t field_id) const;

  const chi_mesh::MeshContinuum& GetSystemCommonGrid() const;

  int64_t MapFieldGlobalIDToSystem(const FieldBlockInfo& info,
                                   int64_t field_global_id) const;
  int64_t MapFieldLocalIDToSystem(const FieldBlockInfo& info,
                                  int64_t field_local_id) const;

  int64_t MapDOF(const chi_mesh::Cell& cell,
                 size_t i,
                 const FieldBlockInfo& info,
                 size_t component_id) const;

  int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                      size_t i,
                      const FieldBlockInfo& info,
                      size_t component_id) const;

  // clang-format off
  class iterator
  {
  public:
    iterator(std::vector<FieldBlockInfo>& ref_vec, size_t ref_element)
      : ref_vec_(ref_vec), ref_element_(ref_element) {}

    iterator operator++() { iterator i = *this; ++ref_element_; return i; }
    iterator operator++(int) { ++ref_element_; return *this; }

    FieldBlockInfo& operator*() { return ref_vec_[ref_element_]; }
    bool operator==(const iterator& rhs) const { return ref_element_ == rhs.ref_element_;;}
    bool operator!=(const iterator& rhs) const { return ref_element_ != rhs.ref_element_; }

  private:
    std::vector<FieldBlockInfo>& ref_vec_;
    size_t ref_element_ = 0;
  };
  // clang-format on

  iterator begin() {return {field_block_info_, 0};}
  iterator end() {return {field_block_info_, field_block_info_.size()};}

  struct FieldBlockInfo
  {
    std::shared_ptr<chi_physics::FieldFunctionGridBased> field_;
    const int64_t num_local_dofs_;
    const int64_t num_global_dofs_;
    const int64_t local_offset_;
    const int64_t ghost_local_offset_;
    const std::vector<int64_t> ghost_ids_;
    const std::vector<std::pair<int64_t, int64_t>> locP_global_block_span_;
    const std::vector<int64_t> locP_system_offsets_;
  };

private:
  static FieldList BuildFieldList(const chi::ParameterBlock& param_array);
  static std::vector<FieldBlockInfo>
  MakeFieldBlockInfo(const FieldList& field_list, int verbosity);
  static int64_t DetermineNumLocalDofs(const std::vector<FieldBlockInfo>&);
  static int64_t DetermineNumGlobalDofs(const std::vector<FieldBlockInfo>&);
  std::vector<int64_t> BuildSystemGhostIDs();

  int verbosity_;
  std::vector<FieldBlockInfo> field_block_info_;
  const int64_t total_num_local_dofs_;
  const int64_t total_num_global_dofs_;
  const std::vector<int64_t> system_ghost_ids_;
};

} // namespace chi_math

#endif // CHITECH_MULTIFIELDCONTAINER_H
