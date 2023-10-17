cmake_minimum_required(VERSION 3.12)

set(TARGET piper)
project(piper LANGUAGES CXX)

#------------------------------------------------ DEPENDENCIES
if (NOT DEFINED CHI_TECH_DIR)
    if (NOT (DEFINED ENV{CHI_TECH_DIR}))
        message(FATAL_ERROR "***** CHI_TECH_DIR is not set *****")
    else()
        set(CHI_TECH_DIR "$ENV{CHI_TECH_DIR}")
    endif()
endif()
message(STATUS "CHI_TECH_DIR set to ${CHI_TECH_DIR}")

include("${CHI_TECH_DIR}/resources/CMakeMacros/Downstream.cmake")

set(SOURCES "main.cc")
add_subdirectory("${PROJECT_SOURCE_DIR}/modules")
add_subdirectory("${PROJECT_SOURCE_DIR}/framework")

include_directories("${PROJECT_SOURCE_DIR}/modules")
include_directories("${PROJECT_SOURCE_DIR}/framework")
include_directories("${PROJECT_SOURCE_DIR}/CoolProp/include")
include_directories("${PROJECT_SOURCE_DIR}/CoolProp/externals/fmtlib/include")

link_directories("${PROJECT_SOURCE_DIR}/CoolProp/build")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")


add_executable(${TARGET} "${SOURCES}")
if(UNIX AND NOT APPLE)
    target_link_libraries(${TARGET} ${CHI_LIBS} -Wl,--whole-archive ChiLib -Wl,--no-whole-archive
            CoolProp)
elseif(APPLE)
    target_link_libraries(${TARGET} ${CHI_LIBS} -Wl,-all_load ChiLib CoolProp)
endif()

file(WRITE ${PROJECT_SOURCE_DIR}/Makefile "subsystem:
" "	$(MAKE) -C chi_build 

"
        "clean:
	$(MAKE) -C chi_build clean
")
