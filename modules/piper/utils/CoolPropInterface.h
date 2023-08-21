#ifndef PIPER_COOLPROPINTERFACE_H
#define PIPER_COOLPROPINTERFACE_H

#include <string>

namespace piper
{

std::string TranslateChiPropertyName2CoolProp(const std::string& chi_name);
double PropSI(const std::string& output_name,
              const std::string& v0n,
              double v0v,
              const std::string& v1n,
              double v1v,
              const std::string& fluid_name);

}

#endif // PIPER_COOLPROPINTERFACE_H
