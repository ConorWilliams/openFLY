cmake_minimum_required(VERSION 3.19)

project(
    OpenKIM-API
    DESCRIPTION "Download, configure, build, and install the OpenKIM API and all its models"
)

# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()


include(ExternalProject)

# -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_FC_COMPILER=${CMAKE_FC_COMPILER}  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    

ExternalProject_Add(
    kim-api
    URL "https://s3.openkim.org/kim-api/kim-api-2.3.0.txz"
    URL_HASH "SHA256=93673bb8fbc0625791f2ee67915d1672793366d10cabc63e373196862c14f991"
    PREFIX ${CMAKE_BINARY_DIR}/kim-api
    CMAKE_ARGS ${CMAKE_BINARY_DIR}/kim-api/src/kim-api  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/kim-api/install
    TEST_COMMAND make test
    TEST_BEFORE_INSTALL TRUE
    STEP_TARGETS install
)

ExternalProject_Add_Step(
    kim-api fetch-models
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/kim-api/install/bin
    COMMAND /bin/bash -c "source kim-api-activate && ./kim-api-collections-management install system EAM_Dynamo_Wen_2021_FeH__MO_634187028437_000"
    COMMENT "Insatalling the KIM models"
)