cmake_minimum_required(VERSION 3.5)

project(DTOKSUproject)

enable_language(CXX)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
#set(CMAKE_CXX_FLAGS -Wall)

option(BUILD_TESTS  "Build test executables" OFF)
option(BUILD_NETCDF  "Build with NetCDF executables" ON)

option(BUILD_DEBUG  "Build with low-level debug" OFF)
option(BUILD_DEEP_DEBUG  "Build with deep-level debug" OFF)
option(BUILD_PAUSE_DEBUG  "Build with PAUSE" OFF)

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/build) 

if(BUILD_NETCDF)
	add_definitions(-DNETCDF_SWITCH)
        
#	set(NETCDF_DIR "/usr/local/netcdf-4.1.3")
#	set(NETCDF_INCLUDES "/usr/local/netcdf-4.1.3/include")
#	set(NETCDF_LIBRARIES "/usr/local/netcdf-4.1.3/lib")
#	set(NETCDF_LIBRARIES_C "/usr/local/netcdf-4.1.3/lib/libnetcdf.so")
#	set(NETCDF_LIBRARIES_CXX "/usr/local/netcdf-4.1.3/lib/libnetcdf_c++.so")
#	set (NETCDF_CXX "YES")
#	set(HDF5_USE_STATIC_LIBRARIES ON)
	find_package (HDF5 REQUIRED COMPONENTS CXX HL)
	find_package (ZLIB REQUIRED)
        find_package (PkgConfig REQUIRED)
        pkg_check_modules (netcdf REQUIRED netcdf)
#        message (${netcdf_DIR})
#        message (${netcdf_INCLUDE_DIRS})
#        message (${netcdf_LIBRARIES})
#        message (${netcdf_LIBRARIES_C})
#        message (${netcdf_LIBRARIES_CXX})
#	find_package (NetCDF REQUIRED)
endif(BUILD_NETCDF)

if(BUILD_DEBUG)
	add_definitions(-DPLASMAGRID_DEBUG)
	add_definitions(-DMODEL_DEBUG)
	add_definitions(-DCHARGING_DEBUG)
	add_definitions(-DHEATING_DEBUG)
	add_definitions(-DFORCE_DEBUG)
	add_definitions(-DMATTER_DEBUG)
	add_definitions(-DELEMENT_DEBUG)
	add_definitions(-DDTOKSU_DEBUG)
	add_definitions(-DDM_DEBUG_DEBUG)
endif(BUILD_DEBUG)

if(BUILD_DEEP_DEBUG)
	add_definitions(-DHEATING_DEEP_DEBUG)
	add_definitions(-DFORCE_DEEP_DEBUG)
	add_definitions(-DMATTER_DEEP_DEBUG)
	add_definitions(-DELEMENT_DEEP_DEBUG)
	add_definitions(-DDTOKSU_DEEP_DEBUG)
endif(BUILD_DEEP_DEBUG)

if(BUILD_PAUSE_DEBUG)
	add_definitions(-DPAUSE)
endif(BUILD_PAUSE_DEBUG)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_executable(dtoksu ${PROJECT_SOURCE_DIR}/main.cpp ${sources} ${headers})

if(BUILD_NETCDF)
	include_directories(${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} ${netcdf_INCLUDE} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/include)
else()
	include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/include)
endif(BUILD_NETCDF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

if(BUILD_TESTS)
	add_subdirectory(Tests)
endif(BUILD_TESTS)

add_library(DTOKSFunc ${PROJECT_SOURCE_DIR}/src/Functions.cpp ${PROJECT_SOURCE_DIR}/src/Constants.cpp ${PROJECT_SOURCE_DIR}/src/threevector.cpp)
add_library(DTOKSCore ${PROJECT_SOURCE_DIR}/src/PlasmaFluxes.cpp ${PROJECT_SOURCE_DIR}/src/CurrentTerms.cpp ${PROJECT_SOURCE_DIR}/src/ForceTerms.cpp ${PROJECT_SOURCE_DIR}/src/HeatTerms.cpp ${PROJECT_SOURCE_DIR}/src/Beryllium.cpp ${PROJECT_SOURCE_DIR}/src/Deuterium.cpp ${PROJECT_SOURCE_DIR}/src/Tungsten.cpp ${PROJECT_SOURCE_DIR}/src/Graphite.cpp ${PROJECT_SOURCE_DIR}/src/Iron.cpp ${PROJECT_SOURCE_DIR}/src/Lithium.cpp ${PROJECT_SOURCE_DIR}/src/Molybdenum.cpp ${PROJECT_SOURCE_DIR}/src/Matter.cpp ${PROJECT_SOURCE_DIR}/src/ChargingModel.cpp ${PROJECT_SOURCE_DIR}/src/HeatingModel.cpp ${PROJECT_SOURCE_DIR}/src/ForceModel.cpp ${PROJECT_SOURCE_DIR}/src/Model.cpp ${PROJECT_SOURCE_DIR}/src/MathHeader.cpp ${PROJECT_SOURCE_DIR}/src/solveMOMLEM.cpp )

if(BUILD_NETCDF)
	target_include_directories(dtoksu PUBLIC ${netcdf_INCLUDE_DIRS} )
	target_link_libraries(dtoksu ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${netcdf_LIBRARIES} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/lib/libconfig4cpp.a)
else()
	target_link_libraries(dtoksu ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/lib/libconfig4cpp.a)
endif(BUILD_NETCDF)

install(TARGETS dtoksu DESTINATION bin)
