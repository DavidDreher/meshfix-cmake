add_library(TMesh edge.cpp  io.cpp  tin.cpp  triangle.cpp  vertex.cpp io_matlab.cpp)
target_include_directories (TMesh PUBLIC ../../include/TMesh)
target_link_libraries(TMesh LINK_PUBLIC Kernel)
target_compile_options(TMesh PRIVATE -fpermissive)

