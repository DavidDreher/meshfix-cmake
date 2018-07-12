add_library(Kernel coordinates.cpp  graph.cpp  heap.cpp  jqsort.cpp  list.cpp  matrix.cpp  orientation.c  point.cpp  tmesh.cpp)
target_include_directories (Kernel PUBLIC ../..//include/Kernel)
target_compile_options(Kernel PRIVATE -fpermissive)

