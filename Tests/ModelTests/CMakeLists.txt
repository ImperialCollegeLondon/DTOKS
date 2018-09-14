file(GLOB testheaders ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/Tests/ModelTests/include/*.hh)

include_directories( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Tests/ModelTests/include)

add_executable (model_test model_test.cpp ${testheaders})

target_link_libraries(model_test DTOKSCore DTOKSFunc )

add_test(NAME MODELTest COMMAND model_test)
