cmake_minimum_required (VERSION 3.1.0)

project(base-scimath VERSION 0.0.1 LANGUAGES C CXX)


if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/askap-cmake)
    execute_process(COMMAND git pull 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/askap-cmake)
else ()	
    execute_process(COMMAND git clone https://bitbucket.csiro.au/scm/askapsdp/askap-cmake.git
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif ()
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/askap-cmake")

option (CXX11 "Compile as C++11 if possible" YES)
option (ENABLE_SHARED "Build shared libraries" YES)
option (ENABLE_RPATH "Include rpath in executables and shared libraries" YES)


# uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()



if (ENABLE_SHARED)
option (BUILD_SHARED_LIBS "" YES)
    if (ENABLE_RPATH)
        # Set RPATH to use for installed targets; append linker search path
        set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
        set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif (ENABLE_RPATH)
endif(ENABLE_SHARED)


# find packages
find_package(Boost REQUIRED COMPONENTS system signals program_options thread chrono )
find_package(LofarCommon REQUIRED)
find_package(log4cxx REQUIRED)
find_package(Casacore REQUIRED COMPONENTS  ms images mirlib coordinates fits lattices measures scimath scimath_f tables casa)
find_package(Askap REQUIRED)
find_package(Imagemath REQUIRED)
find_package(LofarBlob REQUIRED)
find_package(FFTW REQUIRED)
find_package(GSL REQUIRED)
find_package(CPPUnit)

# include directories
include_directories( ${CMAKE_SOURCE_DIR} )
include_directories (${Boost_INCLUDE_DIRS})
add_definitions(-DHAVE_BOOST)
include_directories(${log4cxx_INCLUDE_DIRS})
add_definitions(-DHAVE_LOG4CXX)
include_directories(${LofarCommon_INCLUDE_DIRS})
include_directories(${CASACORE_INCLUDE_DIRS})
add_definitions(-DHAVE_AIPSPP)
include_directories(${ASKAP_INCLUDE_DIRS})
include_directories(${IMAGEMATH_INCLUDE_DIRS})
message(STATUS ${IMAGEMATH_INCLUDE_DIRS})
include_directories(${LofarBlob_INCLUDE_DIRS})
include_directories(${FFTW_INCLUDE_DIRS})

if (CPPUNIT_FOUND)
	include_directories(${CPPUNIT_INCLUDE_DIRS})
endif (CPPUNIT_FOUND)

add_library(askap_scimath

./fitting/ParamsTable.cc
./fitting/ComplexDiffMatrix.cc
./fitting/Quality.cc
./fitting/CompositeEquation.cc
./fitting/Axes.cc
./fitting/NormalEquationsStub.cc
./fitting/PolXProducts.cc
./fitting/ImagingEquation.cc
./fitting/LinearSolver.cc
./fitting/Solveable.cc
./fitting/DesignMatrix.cc
./fitting/GenericEquation.cc
./fitting/Equation.cc
./fitting/PolynomialEquation.cc
./fitting/ParamsCasaTable.cc
./fitting/ImagingNormalEquations.cc
./fitting/Params.cc
./fitting/ComplexDiff.cc
./fitting/ISerializable.cc
./fitting/Solver.cc
./fitting/GenericNormalEquations.cc
./utils/ComplexGaussianNoise.cc
./utils/ImageUtils.cc
./utils/PolConverter.cc
./utils/CasaBlobUtils.cc
./utils/EigenDecompose.cc
./utils/SharedGSLTypes.cc
./utils/SpheroidalFunction.cc
./utils/DelayEstimator.cc
./utils/MultiDimArrayPlaneIter.cc
./utils/PaddingUtils.cc
./utils/MultiDimPosIter.cc
./fft/FFTWrapper.cc
./lsqr_solver/LSQRSolver.cc
./lsqr_solver/MathUtils.cc
./lsqr_solver/ModelDamping.cc
./lsqr_solver/SparseMatrix.cc
)

install (FILES

./fitting/PolynomialEquation.h
./fitting/INormalEquations.h
./fitting/Params.h
./fitting/Equation.h
./fitting/SVDecompose.h
./fitting/NormalEquationsStub.h
./fitting/Solver.h
./fitting/Quality.h
./fitting/GenericNormalEquations.h
./fitting/ComplexDiff.h
./fitting/ImagingEquation.h
./fitting/LinearSolver.h
./fitting/Axes.h
./fitting/ComplexDiffMatrix.h
./fitting/CompositeEquation.h
./fitting/ISerializable.h
./fitting/ParamsTable.h
./fitting/PolXProducts.h
./fitting/GenericEquation.h
./fitting/DesignMatrix.h
./fitting/Solveable.h
./fitting/ParamsCasaTable.h
./fitting/ImagingNormalEquations.h
./fitting/ComplexDiffMatrix.tcc

DESTINATION include/askap/scimath/fitting
)

install (FILES

./utils/SpheroidalFunction.h
./utils/SharedGSLTypes.h
./utils/MultiDimPosIter.h
./utils/CasaBlobUtils.h
./utils/DelayEstimator.h
./utils/ChangeMonitor.h
./utils/ImageUtils.h
./utils/DeepCopyUtils.h
./utils/PaddingUtils.h
./utils/ComplexGaussianNoise.h
./utils/FixedSizeCache.h
./utils/MultiDimArrayPlaneIter.h
./utils/EigenDecompose.h
./utils/PhaseUnwrapper.h
./utils/PolConverter.h
./utils/PaddingUtils.tcc
./utils/PhaseUnwrapper.tcc
./utils/MultiDimArrayPlaneIter.tcc
./utils/CasaBlobUtils.tcc
./utils/FixedSizeCache.tcc

DESTINATION include/askap/scimath/utils
)

install (FILES

./fft/FFTWrapper.h
DESTINATION include/askap/scimath/fft 

)
install (FILES
	
lsqr_solver/GlobalTypedefs.h
lsqr_solver/LSQRSolver.h
lsqr_solver/MathUtils.h
lsqr_solver/ModelDamping.h
lsqr_solver/SparseMatrix.h

DESTINATION include/askap/scimath/lqsr_solver
)

target_link_libraries(
	askap_scimath
	${CASACORE_LIBRARIES}
	${ASKAP_LIBRARY}
	${IMAGEMATH_LIBRARY}
	${log4cxx_LIBRARY}
	${Boost_LIBRARIES}
	${LofarCommon_LIBRARY}
	${LofarBlob_LIBRARY}
	${FFTW_LIBRARIES}
	${GSL_LIBRARIES}
)

install (
TARGETS askap_scimath
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

if (CPPUNIT_FOUND)
	add_executable(
		tutils
		tests/utils/tutils.cc
	)
	target_include_directories(
		tutils
		PRIVATE
		tests/utils/
	)
	target_link_libraries(
		tutils
		askap_scimath
		${CPPUNIT_LIBRARIES}
	)
	install (
	TARGETS tutils
	DESTINATION bin
	)
endif()

