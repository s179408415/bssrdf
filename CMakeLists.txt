cmake_minimum_required(VERSION 3.0.6)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

function(source_group_by_dir source_files)
if(MSVC)
	set(sgbd_cur_dir ${CMAKE_CURRENT_SOURCE_DIR}/src)
	foreach(sgbd_file ${${source_files}})
		string(REGEX REPLACE ${sgbd_cur_dir}/\(.*\) \\1 sgbd_fpath ${sgbd_file})
		string(REGEX REPLACE "\(.*\)/.*" \\1 sgbd_group_name ${sgbd_fpath})
		string(COMPARE EQUAL ${sgbd_fpath} ${sgbd_group_name} sgbd_nogroup)
		string(REPLACE "/" "\\" sgbd_group_name ${sgbd_group_name})
		if(sgbd_nogroup)
			set(sgbd_group_name "\\")
		endif(sgbd_nogroup)
		source_group(${sgbd_group_name} FILES ${sgbd_file})
	endforeach(sgbd_file)
endif(MSVC)
endfunction(source_group_by_dir)
project(bssrdf)
add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake 
	COMMAND conan install .. --build=missing -s build_type=Debug
	DEPENDS ../conanfile.txt)
add_custom_target(conan DEPENDS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake SOURCES ${CMAKE_SOURCE_DIR}/conanfile.txt)
#set_target_properties(t1 PROPERTIES FOLDER Tools)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

file(GLOB_RECURSE INC_FILES src/*.h)
file(GLOB_RECURSE SRC_FILES src/*.cpp)
set(ALL_FILES ${INC_FILES} ${SRC_FILES})
source_group_by_dir(ALL_FILES)


add_executable(bssrdf ${SRC_FILES} ${INC_FILES})
target_link_libraries(bssrdf ${CONAN_LIBS} ${cuda})
add_dependencies(bssrdf conan)
#set_target_properties(bssrdf PROPERTIES FOLDER Engine)
#set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT bssrdf)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT bssrdf)
