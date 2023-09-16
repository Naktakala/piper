#include "math/Systems/EquationSystemTimeData.h"

#include "chi_log.h"

namespace chi_math
{

std::string TimeIDName(TimeID time_id)
{
  switch (time_id)
  {
    case TimeID::T_PLUS_1: return "TIME_T_PLUS_1";
    case TimeID::T: return "TIME_T";
    case TimeID::T_MINUS_1: return "TIME_T_MINUS_1";
    case TimeID::T_MINUS_2: return "TIME_T_MINUS_2";
    case TimeID::T_MINUS_3: return "TIME_T_MINUS_3";
    default:
      ChiLogicalError("Unknown time id");
  }
}

}