cmake_minimum_required(VERSION 3.7)

project(biot_savart)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
#set(CMAKE_SUPPRESS_REGENERATION true)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

# Find OpenMP
find_package(OpenMP)

############################################################
# Create a library
############################################################

#Generate the shared library from the library sources
add_library(biot_savart_shared SHARED biot_savart.cpp) 
set_target_properties(biot_savart_shared PROPERTIES OUTPUT_NAME biot_savart)

if (UNIX)
    target_link_libraries(biot_savart_shared pthread)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64" )
    install(TARGETS biot_savart_shared LIBRARY DESTINATION lib)    
endif()  

if (WIN32) 
    add_compile_definitions(__EXPORT_CLASS_BIOTSAVART__) 
endif() 

find_package(Matlab)
if(Matlab_FOUND)
 	matlab_add_mex(NAME biot_savart SRC biot_savart.cpp R2018a) 
endif()


