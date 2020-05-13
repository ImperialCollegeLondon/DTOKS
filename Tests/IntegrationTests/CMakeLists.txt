file(GLOB testheaders ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/Tests/IntegrationTests/include/*.hh)

include_directories( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Tests/IntegrationTests/include ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/include)

add_executable (int_test integration_test.cpp ${testheaders})

target_link_libraries(int_test DTOKSCore DTOKSFunc ${PROJECT_SOURCE_DIR}/Dependencies/config4cpp/lib/libconfig4cpp.a)

add_test(NAME INTTest COMMAND int_test)
