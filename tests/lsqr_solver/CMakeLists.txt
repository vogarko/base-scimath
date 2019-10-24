add_executable(tlsqr_solver tlsqr_solver.cc)
target_link_libraries(tlsqr_solver
	askap_scimath
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tlsqr_solver
	COMMAND tlsqr_solver
	)
