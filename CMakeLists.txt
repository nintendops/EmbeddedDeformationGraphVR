project(DeformationVR LANGUAGES CXX)

cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# If not set, determines the running platform architecture.
if(NOT PLATFORM)
  if(CMAKE_SIZEOF_VOID_P MATCHES 8)
    set(PLATFORM 64)
  else()
    set(PLATFORM 32)
  endif()
endif()
message(STATUS "Compilation set for ${PLATFORM}bits architectures.")

set(SDL_REQUIRED_LIBRARIES SDL2)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions(-DLINUX -DPOSIX)
  set(ARCH_TARGET linux64)

  if(${PLATFORM} MATCHES 32)
    message(WARNING "OpenVR x86 binaries not provided on GNU/Linux.")
  endif()
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_MACOSX_RPATH 0)
  add_definitions(-DOSX -DPOSIX)
  set(ARCH_TARGET osx32)

elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(SDL_REQUIRED_LIBRARIES ${SDL_REQUIRED_LIBRARIES} SDL2main)
  add_definitions(-D_WIN32)
  set(ARCH_TARGET win${PLATFORM})

  # Binaries path for thirdparties are not generics so we try to guess their suffixes.
  set(WINDOWS_PATH_SUFFIXES win${PLATFORM} Win${PLATFORM} x${PLATFORM})

  if(${PLATFORM} MATCHES 64)
    message(WARNING "SDL x64 runtime binaries not provided on Windows.")
  endif()
endif()

# -----------------------------------------------------------------------------
## PATHS ##

set(THIRDPARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
set(SHARED_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/shared)
set(GMM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gmm)

# Check that the steamVR SDK is installed 
# (needed to prevent a segfault in OpenVR).
if(CMAKE_HOST_UNIX)
  find_file(OPENVRPATHS openvrpaths.vrpath PATHS $ENV{HOME}/.config/openvr "$ENV{HOME}/Library/Application Support/OpenVR/.openvr")
  if(${OPENVRPATHS} MATCHES OPENVRPATHS-NOTFOUND)
    message(FATAL_ERROR "${OPENVRPATHS} Please install SteamVR SDK to continue..")
  endif()
endif()

# setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# setup GLUT
find_package(GLUT)
include_directories(${GLUT_INCLUDE_DIR})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found")
endif()

# setup OpenGL
find_package(OpenGL)
include_directories(${OpenGL_INCLUDE_DIRS})
if(NOT OPENGL_FOUND)
    message(ERROR " OpenGL not found")
endif()

# setup OpenMesh
find_package(OpenMesh)
include_directories(${OPENMESH_INCLUDE_DIRS})
if(NOT OPENMESH_FOUND)
    message(ERROR " OpenMesh not found")
endif()

# setup Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# setup OpenVR API
find_library(OPENVR_LIBRARIES
  NAMES
    openvr_api
  PATHS
    ${CMAKE_CURRENT_SOURCE_DIR}/OpenVR/bin
    ${CMAKE_CURRENT_SOURCE_DIR}/OpenVR/lib
  PATH_SUFFIXES
    osx32
    linux64
    ${WINDOWS_PATH_SUFFIXES}
  NO_DEFAULT_PATH
  NO_CMAKE_FIND_ROOT_PATH
)
set(OPENVR_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/OpenVR/headers)


## OpenGL / GLU
find_package(OpenGL REQUIRED)

## GLEW 1.11
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # On GNU/Linux, glew can be found on the package manager. 
  find_package(GLEW 1.11 REQUIRED)
else()
  # Otherwise, use the shipped binaries.
  find_library(GLEW_LIBRARIES
    NAMES
      GLEW
      glew32
    PATHS
      ${THIRDPARTY_DIR}/glew/glew-1.11.0/lib/Release
    PATH_SUFFIXES
      osx32
      ${WINDOWS_PATH_SUFFIXES}
  )
  set(GLEW_INCLUDE_DIR ${THIRDPARTY_DIR}/glew/glew-1.11.0/include)
endif()

## SDL 2
## osx32 and linux32 cmake can spot each others' .so and .a files, so keep them apart
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  foreach(lib ${SDL_REQUIRED_LIBRARIES})
    find_library(${lib}_LIBRARY
      NAMES
        ${lib}
      PATHS
        ${THIRDPARTY_DIR}/sdl2-2.0.3/bin
      PATH_SUFFIXES
        osx32
    )
    list(APPEND SDL2_LIBRARIES ${${lib}_LIBRARY})
  endforeach()
else()
  foreach(lib ${SDL_REQUIRED_LIBRARIES})
    find_library(${lib}_LIBRARY
      NAMES
        ${lib}
      PATHS
        ${THIRDPARTY_DIR}/sdl2-2.0.3/bin
      PATH_SUFFIXES
        linux${PLATFORM}
        ${WINDOWS_PATH_SUFFIXES}
    )
    list(APPEND SDL2_LIBRARIES ${${lib}_LIBRARY})
  endforeach()
endif()
set(SDL2_INCLUDE_DIR ${THIRDPARTY_DIR}/sdl2-2.0.3/include)


# required by OpenMesh
set_property(
    DIRECTORY
    APPEND PROPERTY COMPILE_DEFINITIONS _USE_MATH_DEFINES
)

# collect sources
file(GLOB deformer_sources deformation/*.cc)
file(GLOB deformer_headers deformation/*.hh)

file(GLOB viewer_sources viewer/*.cc)
file(GLOB viewer_headers viewer/*.hh)

file(GLOB ann_sources ANN/*.cpp)
file(GLOB ann_headers ANN/*.h)

file(GLOB_RECURSE isoex_sources isoex/*.cc)
file(GLOB_RECURSE isoex_headers isoex/*.hh)

file(GLOB SHARED_SRC_FILES
  ${SHARED_SRC_DIR}/*.cpp
  ${SHARED_SRC_DIR}/*.h
)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${OPENVR_INCLUDE_DIR}
	${GMM_DIR}
	${OPENGL_INCLUDE_DIR}
	${GLEW_INCLUDE_DIR}
	${SDL2_INCLUDE_DIR}
	${OPENVR_INCLUDE_DIR}
)

add_definitions(-DINCLUDE_TEMPLATES)

# compile and link
add_library(ann ${ann_sources} ${ann_headers})
add_library(isoex ${isoex_sources} ${isoex_headers})
add_executable(edVR ${deformer_sources} ${SHARED_SRC_FILES} ${deformer_headers})
add_executable(viewer ${viewer_sources} ${viewer_headers})
target_link_libraries(viewer ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OPENMESH_LIBRARIES})
target_link_libraries(edVR 
	${OPENMESH_LIBRARIES} 
	${OPENGL_LIBRARIES} 
	${OPENVR_LIBRARIES} 
	${GLEW_LIBRARIES}
	${SDL2_LIBRARIES}
	${CMAKE_DL_LIBS}
	isoex ann Eigen3::Eigen)

target_compile_features(edVR PUBLIC cxx_std_11)
#CUDA requirements
#set_target_properties(edVR PROPERTIES CUDA_SEPARABLE_COMPILATION ON)