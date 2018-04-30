cmake_minimum_required(VERSION 3.5)

project(DTOKSUproject)

enable_language(CXX)

set (CMAKE_CXX_STANDARD 14)

set (NETCDF_CXX "YES")
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/build) 

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(main ${PROJECT_SOURCE_DIR}/main.cpp ${sources} ${headers})

find_package (HDF5 REQUIRED COMPONENTS CXX HL)
find_package (ZLIB REQUIRED)
find_package (NetCDF REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} ${NETCDF_INCLUDES} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/include)


target_link_libraries(main ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${NETCDF_LIBRARIES_CXX} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/lib/libconfig4cpp.a)

install(TARGETS main DESTINATION bin)
