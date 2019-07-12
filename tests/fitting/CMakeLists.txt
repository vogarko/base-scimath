add_executable(tfitting tfitting.cc)
target_link_libraries(tfitting
	askap_scimath
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tfitting
	COMMAND tfitting
	)
