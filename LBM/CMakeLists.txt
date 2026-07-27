cmake_minimum_required(VERSION 3.20)


project(
    ISLBM_Cylinder
    VERSION 1.0
    LANGUAGES CXX
)



set(
    CMAKE_CXX_STANDARD
    20
)



find_package(OpenMP REQUIRED)



add_executable(

    islbm


    app/main.cpp


    src/Mesh.cpp
    src/Geometry.cpp
    src/Boundary.cpp
    src/Collision.cpp
    src/Streaming.cpp
    src/LBM.cpp
    src/VTKWriter.cpp

)



target_include_directories(

    islbm

    PRIVATE

    include

)



target_link_libraries(

    islbm

    PRIVATE

    OpenMP::OpenMP_CXX

)