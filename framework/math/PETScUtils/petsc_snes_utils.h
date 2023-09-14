#ifndef CHI_PETSC_SNES_UTILS_H
#define CHI_PETSC_SNES_UTILS_H

#include <petscsnes.h>

namespace chi_math::PETScUtils
{

PetscErrorCode
BasicSNESMonitor(SNES snes, PetscInt iter, PetscReal rnorm, void*);

}

#endif // CHI_PETSC_SNES_UTILS_H
