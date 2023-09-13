#ifndef PIPER_PETSC_SNES_UTILS_H
#define PIPER_PETSC_SNES_UTILS_H

#include <petscsnes.h>

namespace chi_math::PETScUtils
{

PetscErrorCode
BasicSNESMonitor(SNES snes, PetscInt iter, PetscReal rnorm, void*);

}

#endif // PIPER_PETSC_SNES_UTILS_H
