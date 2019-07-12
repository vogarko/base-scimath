add_executable(tfft tfft.cc)
target_link_libraries(tfft
	askap_scimath
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tfft
	COMMAND tfft
	)
