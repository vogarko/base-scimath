add_executable(tutils tutils.cc)
target_link_libraries(tutils
	askap_scimath
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tutils
	COMMAND tutils
	)
