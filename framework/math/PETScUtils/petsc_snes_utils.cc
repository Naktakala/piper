#include "petsc_snes_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>

namespace chi_math::PETScUtils
{

PetscErrorCode
BasicSNESMonitor(SNES snes, PetscInt iter, PetscReal rnorm, void*)
{
  const char* snes_name;
  SNESGetOptionsPrefix(snes, &snes_name);

  //Default to this if ksp_name is NULL
  const char NONAME_SOLVER[] = "NoName-Solver\0";

  if (snes_name == nullptr)
    snes_name = NONAME_SOLVER;

  //Print message
  std::stringstream buff;
  buff
    << snes_name
    << " non-linear iteration "
    << std::setw(4) << iter
    << " - Residual "
    << std::scientific << std::setprecision(7) << rnorm
    << std::endl;

  Chi::log.Log() << buff.str();

  return 0;
}

}