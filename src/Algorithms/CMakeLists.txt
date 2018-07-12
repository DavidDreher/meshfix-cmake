add_library (Algorithms checkAndRepair.cpp  detectIntersections.cpp  holeFilling.cpp  marchIntersections.cpp  
subdivision.cpp)
target_link_libraries (Algorithms LINK_PUBLIC TMesh)
target_compile_options(Algorithms PRIVATE -fpermissive)

