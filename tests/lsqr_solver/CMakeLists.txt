add_executable(tlsqr_solver tlsqr_solver.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tlsqr_solver
	askap::scimath
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tlsqr_solver
	COMMAND tlsqr_solver
	)
