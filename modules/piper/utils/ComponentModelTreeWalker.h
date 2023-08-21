#ifndef PIPER_COMPONENTMODELTREEWALKER_H
#define PIPER_COMPONENTMODELTREEWALKER_H

#include <map>
#include <string>
#include <memory>
#include <functional>
#include <set>

namespace piper
{

class ComponentModel;

class ComponentModelTreeWalker
{
public:
  typedef std::vector<std::unique_ptr<ComponentModel>> ModelMap;
  typedef std::function<void(ComponentModel&)> OperatorFunction;

  ComponentModelTreeWalker(ModelMap& model_map,
                           size_t starting_component_id,
                           const OperatorFunction& operator_function);

  void Execute();

  static void Execute(ModelMap& model_map,
                      size_t starting_component_id,
                      const OperatorFunction& operator_function);

protected:
  void Walk(ComponentModel& component, std::set<size_t>& comp_ids_visited);

  ModelMap& model_map_;
  const size_t starting_component_id_;
  const OperatorFunction& op_;
};

} // namespace piper

#endif // PIPER_COMPONENTMODELTREEWALKER_H
