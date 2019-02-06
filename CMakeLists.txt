cmake_minimum_required(VERSION 3.5)

project(DTOKSUproject)

enable_language(CXX)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

option(BUILD_TESTS  "Build test executables" OFF)
option(BUILD_NETCDF  "Build with NetCDF executables" ON)

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/build) 

if(BUILD_NETCDF)
	add_definitions(-DNETCDF_SWITCH)
	set (NETCDF_CXX "YES")
	find_package (HDF5 REQUIRED COMPONENTS CXX HL)
	find_package (ZLIB REQUIRED)
	find_package (NetCDF REQUIRED)
endif(BUILD_NETCDF)



file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(dtoksu ${PROJECT_SOURCE_DIR}/main.cpp ${sources} ${headers})

if(BUILD_NETCDF)
	include_directories(${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} ${NETCDF_INCLUDES} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/include)
else()
	include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/include)
endif(BUILD_NETCDF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

if(BUILD_TESTS)
	add_subdirectory(Tests)
endif(BUILD_TESTS)

add_library(DTOKSFunc ${PROJECT_SOURCE_DIR}/src/Functions.cpp ${PROJECT_SOURCE_DIR}/src/Constants.cpp ${PROJECT_SOURCE_DIR}/src/threevector.cpp)
add_library(DTOKSCore ${PROJECT_SOURCE_DIR}/src/Beryllium.cpp ${PROJECT_SOURCE_DIR}/src/Deuterium.cpp ${PROJECT_SOURCE_DIR}/src/Tungsten.cpp ${PROJECT_SOURCE_DIR}/src/Graphite.cpp ${PROJECT_SOURCE_DIR}/src/Iron.cpp ${PROJECT_SOURCE_DIR}/src/Lithium.cpp ${PROJECT_SOURCE_DIR}/src/Molybdenum.cpp ${PROJECT_SOURCE_DIR}/src/Matter.cpp ${PROJECT_SOURCE_DIR}/src/ChargingModel.cpp ${PROJECT_SOURCE_DIR}/src/HeatingModel.cpp ${PROJECT_SOURCE_DIR}/src/ForceModel.cpp ${PROJECT_SOURCE_DIR}/src/Model.cpp ${PROJECT_SOURCE_DIR}/src/MathHeader.cpp ${PROJECT_SOURCE_DIR}/src/solveMOMLEM.cpp )

if(BUILD_NETCDF)
	target_link_libraries(dtoksu ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${NETCDF_LIBRARIES_CXX} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/lib/libconfig4cpp.a)
else()
	target_link_libraries(dtoksu ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/lib/libconfig4cpp.a)
endif(BUILD_NETCDF)

install(TARGETS dtoksu DESTINATION bin)
