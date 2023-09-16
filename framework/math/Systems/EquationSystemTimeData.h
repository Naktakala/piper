#ifndef CHITECH_EQUATIONSYSTEMTIMEDATA_H
#define CHITECH_EQUATIONSYSTEMTIMEDATA_H

#include <string>

namespace chi_math
{

enum class TimeID : int
{
  T_PLUS_1 = -1,
  T = 0,
  T_MINUS_1 = 1,
  T_MINUS_2 = 2,
  T_MINUS_3 = 3
};

std::string TimeIDName(TimeID time_id);

/**Structure for passing time data.*/
struct EquationSystemTimeData
{
  double dt_ = 0.001;
  double time_ = 0.0;
};

} // namespace chi_math

#endif // CHITECH_EQUATIONSYSTEMTIMEDATA_H
