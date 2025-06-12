cmake_minimum_required(VERSION 3.2 FATAL_ERROR)  

project(muTruck)

find_package(Geant4 REQUIRED ui_all vis_all) 

include(${Geant4_USE_FILE}) 
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc) 
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB MACRO_FILES "*.mac")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

file(GLOB DATA_FILES "*.dat")
file(COPY ${DATA_FILES} DESTINATION ${PROJECT_BINARY_DIR})

add_executable(muonTruck muonTruck.cc ${sources} ${headers}) 
target_include_directories(muonTruck PRIVATE include)
target_link_libraries(muonTruck ${Geant4_LIBRARIES})  

add_executable(muonTruck_multithread muonTruck_multithread.cc ${sources} ${headers}) 
target_include_directories(muonTruck_multithread PRIVATE include)
target_link_libraries(muonTruck_multithread ${Geant4_LIBRARIES})  

add_custom_target(muTruck DEPENDS muonTruck)  
