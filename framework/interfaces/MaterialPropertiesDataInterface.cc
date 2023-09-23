#include "MaterialPropertiesDataInterface.h"

#include "materials/MaterialPropertiesData.h"

namespace chi
{

InputParameters MaterialPropertiesDataInterface::GetInputParameters()
{
  InputParameters params;

  params.AddOptionalParameter("material_properties",
                              0,
                              "Handle to a MaterialPropertiesData-object "
                              "containing the material properties"
                              " relevant to the simulation.");

  return params;
}

MaterialPropertiesDataInterface::MaterialPropertiesDataInterface(
  const InputParameters& params)
  : material_properties_data_(GetDataBlock(params))
{
}

std::shared_ptr<MaterialPropertiesData>
MaterialPropertiesDataInterface::GetDataBlock(const InputParameters& params)
{
  if (params.ParametersAtAssignment().Has("material_properties"))
    return Chi::GetStackItemPtrAsType<chi::MaterialPropertiesData>(
      /*stack=*/Chi::object_stack,
      /*handle=*/params.GetParamValue<size_t>("material_properties"),
      /*calling_function_name=*/__FUNCTION__);
  else
    return MaterialPropertiesData::MakeEmptyData();
}

} // namespace chi