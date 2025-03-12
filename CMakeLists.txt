cmake_minimum_required(VERSION 3.10)
project(TimoshenkoEig VERSION 0.1 DESCRIPTION forvscode)
# project  name = executable name

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED True)


option(USE_TIMOSHENKO_EIG "User experiment" ON)

configure_file(NewbieConfig.h.in NewbieConfig.h)

# subdirectory name = library name
if(USE_TIMOSHENKO_EIG)
  add_subdirectory(src)
  list(APPEND EXTRA_LIBS src)
  list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/src")
endif()

add_executable(TimoshenkoEig arclengthfig.cxx)
#set_property(TARGET PROPERTY CXX_STANDARD 14)
#add_executable(circular circular_cross_section.cxx)
#add_executable(rectangular rectangular_cross_section.cxx)


target_link_libraries(TimoshenkoEig PUBLIC ${EXTRA_LIBS})

# add the binary tree to the search path for include files
# so that we will find NewbieConfig.h
target_include_directories(TimoshenkoEig PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

