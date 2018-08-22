file(GLOB testheaders ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/Tests/UnitTests/include/*.hh)

include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Tests/UnitTests/include)

add_executable (unit_test unit_test.cpp ${testheaders})

target_link_libraries(unit_test DTOKSFunc DTOKSCore)

add_test(NAME UNITTest COMMAND unit_test)
