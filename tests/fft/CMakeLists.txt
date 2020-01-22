add_executable(tfft tfft.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tfft
	askap::scimath
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tfft
	COMMAND tfft
	)
