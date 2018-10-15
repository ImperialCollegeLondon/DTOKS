cmake_minimum_required(VERSION 3.5)

project(DTOKSUproject)

enable_language(CXX)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set (NETCDF_CXX "YES")
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/build) 

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(dtoksu ${PROJECT_SOURCE_DIR}/main.cpp ${sources} ${headers})

find_package (HDF5 REQUIRED COMPONENTS CXX HL)
find_package (ZLIB REQUIRED)
find_package (NetCDF REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} ${NETCDF_INCLUDES} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/include)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_subdirectory(Tests)

add_library(DTOKSFunc ${PROJECT_SOURCE_DIR}/src/Functions.cpp ${PROJECT_SOURCE_DIR}/src/Constants.cpp ${PROJECT_SOURCE_DIR}/src/threevector.cpp)
add_library(DTOKSCore ${PROJECT_SOURCE_DIR}/src/Beryllium.cpp ${PROJECT_SOURCE_DIR}/src/Deuterium.cpp ${PROJECT_SOURCE_DIR}/src/Tungsten.cpp ${PROJECT_SOURCE_DIR}/src/Graphite.cpp ${PROJECT_SOURCE_DIR}/src/Iron.cpp ${PROJECT_SOURCE_DIR}/src/Matter.cpp ${PROJECT_SOURCE_DIR}/src/ChargingModel.cpp ${PROJECT_SOURCE_DIR}/src/HeatingModel.cpp ${PROJECT_SOURCE_DIR}/src/ForceModel.cpp ${PROJECT_SOURCE_DIR}/src/Model.cpp ${PROJECT_SOURCE_DIR}/src/MathHeader.cpp)

target_link_libraries(dtoksu ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${NETCDF_LIBRARIES_CXX} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/lib/libconfig4cpp.a)

install(TARGETS dtoksu DESTINATION bin)