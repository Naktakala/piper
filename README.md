# piper
ChiTech based piping simulator

## How to build CoolProp
git submodule update --init CoolProp
cd CoolProp
mkdir build
cd build
cmake .. -DCOOLPROP_STATIC_LIBRARY=ON
cmake --build .
