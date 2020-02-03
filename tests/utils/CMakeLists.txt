add_executable(tutils tutils.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tutils
	askap::scimath
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tutils
	COMMAND tutils
	)
