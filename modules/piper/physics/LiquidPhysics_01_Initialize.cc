#include "LiquidPhysics.h"

#include "piper/MeshGenerators/PiperMeshGenerator.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace piper
{

void LiquidPhysics::Initialize()
{
  Chi::log.Log() << "Making mesh";
  Chi::mpi.Barrier();
  this->MakeMesh();
  Chi::log.Log() << "Initializing Unknowns";
  Chi::mpi.Barrier();
  this->InitializeUnknowns();
}

}