add_executable(tfitting tfitting.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tfitting
	askap::scimath
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tfitting
	COMMAND tfitting
	)
