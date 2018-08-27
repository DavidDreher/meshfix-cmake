find_package(Matlab REQUIRED)
include_directories(${Matlab_INCLUDE_DIRS})

set(CMAKE_SHARED_LIBRARY_SUFFIX ".${Matlab_MEX_EXTENSION}")
set(CMAKE_SHARED_LIBRARY_PREFIX )
set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)

set(LIBS)
list(APPEND LIBS ${Matlab_LIBRARIES})

add_library(mat_meshfix SHARED mat_meshfix.cpp)
target_link_libraries(mat_meshfix LINK_PUBLIC Algorithms TMesh ${LIBS})

target_compile_options(mat_meshfix PRIVATE -fpermissive)


add_custom_command(TARGET mat_meshfix POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mat_meshfix> ${PROJECT_SOURCE_DIR}
        )



#
#add_executable (meshfix meshfix.cpp)
#target_link_libraries (meshfix LINK_PUBLIC Algorithms TMesh)
#
#
#add_library(TMesh edge.cpp  io.cpp  tin.cpp  triangle.cpp  vertex.cpp io_matlab.cpp)
#target_include_directories (TMesh PUBLIC ../../include/TMesh)
#target_link_libraries(TMesh LINK_PUBLIC Kernel)
#target_compile_options(TMesh PRIVATE -fpermissive)