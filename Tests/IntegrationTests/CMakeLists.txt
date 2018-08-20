enable_language(CXX)

set (CMAKE_CXX_STANDARD 14)
set (NETCDF_CXX "YES")
find_package (HDF5 REQUIRED COMPONENTS CXX HL)
find_package (ZLIB REQUIRED)
find_package (NetCDF REQUIRED)

file(GLOB testheaders ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/Tests/IntegrationTests/include/*.hh)

include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Tests/IntegrationTests/include ${HDF5_INCLUDE_DIRS} ${NETCDF_INCLUDES} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/include)

add_executable (int_test integration_test.cpp ${testheaders})

target_link_libraries(int_test DTOKSCore DTOKSFunc ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${NETCDF_LIBRARIES_CXX} ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/lib/libconfig4cpp.a)

add_test(NAME INTTest COMMAND int_test)
