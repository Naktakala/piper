#ifndef PIPER_CONNECTIONINTERFACE_H
#define PIPER_CONNECTIONINTERFACE_H

#include <vector>
#include <cstddef>

namespace piper
{

class ComponentModel;

namespace utils
{
class Connection;
}

class ConnectionInterface
{
public:
  explicit ConnectionInterface(ComponentModel& model);

  struct Info
  {
    size_t connection_index_;
    ComponentModel& component_model_;
    const utils::Connection& connection_;
  };

  // clang-format off
  /**This iterator allows us to iterate over the connections whilst
  * directly getting the components connected to it.*/
  class iterator
  {
  public:
    iterator(ComponentModel& model, size_t i);

    iterator operator++();
    iterator operator++(int);

    Info operator*();
    bool operator==(const iterator& rhs) const;
    bool operator!=(const iterator& rhs) const;

  private:
    ComponentModel& model_;
    const std::vector<utils::Connection>& connections_;
    size_t ref_elem_;
  };
  // clang-format on

  iterator begin();
  iterator end();

  size_t size();

private:
  ComponentModel& component_model_;
};

}

#endif // PIPER_CONNECTIONINTERFACE_H
