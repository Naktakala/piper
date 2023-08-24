#include "ConnectionInterface.h"

#include "piper/utils/Connections.h"
#include "piper/models/ComponentModel.h"

namespace piper
{

ConnectionInterface::ConnectionInterface(ComponentModel& model)
  : component_model_(model)
{
}

ConnectionInterface::iterator::iterator(ComponentModel& model, size_t i)
  : model_(model), connections_(model.ConnectionPoints()), ref_elem_(i)
{
}

ConnectionInterface::iterator ConnectionInterface::iterator::operator++()
{
  iterator i = *this;
  ++ref_elem_;
  return i;
}

ConnectionInterface::iterator ConnectionInterface::iterator::operator++(int)
{
  ++ref_elem_;
  return *this;
}

ConnectionInterface::Info
ConnectionInterface::iterator::operator*()
{
  const auto& connections = model_.ConnectionPoints();
  const auto& connection = connections.at(ref_elem_);
  auto& ref_model = model_.Family().at(connection.connected_comp_id_);
  return Info{ref_elem_, *ref_model, model_.ConnectionPoints().at(ref_elem_)};
}

bool ConnectionInterface::iterator::operator==(const iterator& rhs) const
{
  return ref_elem_ == rhs.ref_elem_;
}

bool ConnectionInterface::iterator::operator!=(const iterator& rhs) const
{
  return ref_elem_ != rhs.ref_elem_;
}

ConnectionInterface::iterator ConnectionInterface::begin()
{
  return iterator(component_model_, 0);
}

ConnectionInterface::iterator ConnectionInterface::end()
{
  auto& connections = component_model_.ConnectionPoints();
  return iterator(component_model_, connections.size());
}

size_t ConnectionInterface::size()
{
  return component_model_.ConnectionPoints().size();
}

} // namespace piper