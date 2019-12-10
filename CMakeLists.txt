cmake_minimum_required (VERSION 3.1.0)

project(base-scimath VERSION 0.0.1 LANGUAGES C CXX)


if (EXISTS ${CMAKE_SOURCE_DIR}/askap-cmake)
    execute_process(COMMAND git pull 
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/askap-cmake)
else ()	
    execute_process(COMMAND git clone https://bitbucket.csiro.au/scm/askapsdp/askap-cmake.git
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif ()
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/askap-cmake")

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

# Yanda packages
include(install_yanda_depends)

set (EXTERNAL_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
find_package(LofarCommon)
if (NOT LofarCommon_FOUND)
	install_yanda_depends(lofar-common)
	find_package(LofarCommon REQUIRED)
endif ()

find_package(LofarBlob)
if (NOT LofarBlob_FOUND)
	install_yanda_depends(lofar-blob)
	find_package(LofarBlob REQUIRED)
endif ()


find_package(Askap)
if (NOT ASKAP_FOUND)
	install_yanda_depends(base-askap)
	find_package(Askap REQUIRED)
endif ()

find_package(Imagemath)
if (NOT IMAGEMATH_FOUND)
	install_yanda_depends(base-imagemath)
	find_package(Imagemath REQUIRED)
endif ()


find_package(AskapParallel)
if (NOT ASKAPPARALLEL_FOUND)
	install_yanda_depends(base-askapparallel)
	find_package(AskapParallel REQUIRED)
endif ()


# find packages
find_package(Boost REQUIRED COMPONENTS system program_options thread chrono )
find_package(log4cxx REQUIRED)
find_package(Casacore REQUIRED COMPONENTS  ms images mirlib coordinates fits lattices measures scimath scimath_f tables casa)
find_package(FFTW REQUIRED)
find_package(GSL REQUIRED)
find_package(CPPUnit)
find_package(MPI REQUIRED)

if (MPI_FOUND)
    add_definitions(-DHAVE_MPI)
    include_directories(${MPI_INCLUDE_PATH})
endif()


# include directories
include_directories (${Boost_INCLUDE_DIRS})
add_definitions(-DHAVE_BOOST)
include_directories(${log4cxx_INCLUDE_DIRS})
add_definitions(-DHAVE_LOG4CXX)
include_directories(${CASACORE_INCLUDE_DIRS})
add_definitions(-DHAVE_AIPSPP)
add_definitions(-Dcasa=casacore)
include_directories(${FFTW_INCLUDES})


if (NOT GSL_VERSION VERSION_LESS 2.0)
	add_definitions(-DHAVE_GSL2)
endif()

if (CASACORE3 OR CXX11)
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if (CPPUNIT_FOUND)
	include_directories(${CPPUNIT_INCLUDE_DIRS})
endif (CPPUNIT_FOUND)

add_library(askap_scimath

./askap/scimath/fitting/ParamsTable.cc
./askap/scimath/fitting/ComplexDiffMatrix.cc
./askap/scimath/fitting/Quality.cc
./askap/scimath/fitting/CompositeEquation.cc
./askap/scimath/fitting/Axes.cc
./askap/scimath/fitting/NormalEquationsStub.cc
./askap/scimath/fitting/PolXProducts.cc
./askap/scimath/fitting/ImagingEquation.cc
./askap/scimath/fitting/LinearSolver.cc
./askap/scimath/fitting/LinearSolverLsqrUtils.cc
./askap/scimath/fitting/LinearSolverUtils.cc
./askap/scimath/fitting/Solveable.cc
./askap/scimath/fitting/DesignMatrix.cc
./askap/scimath/fitting/GenericEquation.cc
./askap/scimath/fitting/Equation.cc
./askap/scimath/fitting/PolynomialEquation.cc
./askap/scimath/fitting/ParamsCasaTable.cc
./askap/scimath/fitting/ImagingNormalEquations.cc
./askap/scimath/fitting/Params.cc
./askap/scimath/fitting/ComplexDiff.cc
./askap/scimath/fitting/ISerializable.cc
./askap/scimath/fitting/Solver.cc
./askap/scimath/fitting/GenericNormalEquations.cc
./askap/scimath/utils/ComplexGaussianNoise.cc
./askap/scimath/utils/ImageUtils.cc
./askap/scimath/utils/PolConverter.cc
./askap/scimath/utils/CasaBlobUtils.cc
./askap/scimath/utils/EigenDecompose.cc
./askap/scimath/utils/SharedGSLTypes.cc
./askap/scimath/utils/SpheroidalFunction.cc
./askap/scimath/utils/DelayEstimator.cc
./askap/scimath/utils/MultiDimArrayPlaneIter.cc
./askap/scimath/utils/PaddingUtils.cc
./askap/scimath/utils/MultiDimPosIter.cc
./askap/scimath/fft/FFTWrapper.cc
./askap/scimath/lsqr_solver/LSQRSolver.cc
./askap/scimath/lsqr_solver/MathUtils.cc
./askap/scimath/lsqr_solver/ModelDamping.cc
./askap/scimath/lsqr_solver/ParallelTools.cc
./askap/scimath/lsqr_solver/SparseMatrix.cc
)

install (FILES

./askap/scimath/fitting/PolynomialEquation.h
./askap/scimath/fitting/INormalEquations.h
./askap/scimath/fitting/Params.h
./askap/scimath/fitting/Equation.h
./askap/scimath/fitting/SVDecompose.h
./askap/scimath/fitting/NormalEquationsStub.h
./askap/scimath/fitting/Solver.h
./askap/scimath/fitting/Quality.h
./askap/scimath/fitting/GenericNormalEquations.h
./askap/scimath/fitting/ComplexDiff.h
./askap/scimath/fitting/ImagingEquation.h
./askap/scimath/fitting/LinearSolver.h
./askap/scimath/fitting/LinearSolverLsqrUtils.cc
./askap/scimath/fitting/LinearSolverUtils.cc
./askap/scimath/fitting/Axes.h
./askap/scimath/fitting/ComplexDiffMatrix.h
./askap/scimath/fitting/CompositeEquation.h
./askap/scimath/fitting/ISerializable.h
./askap/scimath/fitting/ParamsTable.h
./askap/scimath/fitting/PolXProducts.h
./askap/scimath/fitting/GenericEquation.h
./askap/scimath/fitting/DesignMatrix.h
./askap/scimath/fitting/Solveable.h
./askap/scimath/fitting/ParamsCasaTable.h
./askap/scimath/fitting/ImagingNormalEquations.h
./askap/scimath/fitting/ComplexDiffMatrix.tcc

DESTINATION include/askap/scimath/fitting
)

install (FILES

./askap/scimath/utils/SpheroidalFunction.h
./askap/scimath/utils/SharedGSLTypes.h
./askap/scimath/utils/MultiDimPosIter.h
./askap/scimath/utils/CasaBlobUtils.h
./askap/scimath/utils/DelayEstimator.h
./askap/scimath/utils/ChangeMonitor.h
./askap/scimath/utils/ImageUtils.h
./askap/scimath/utils/DeepCopyUtils.h
./askap/scimath/utils/PaddingUtils.h
./askap/scimath/utils/ComplexGaussianNoise.h
./askap/scimath/utils/FixedSizeCache.h
./askap/scimath/utils/MultiDimArrayPlaneIter.h
./askap/scimath/utils/EigenDecompose.h
./askap/scimath/utils/PhaseUnwrapper.h
./askap/scimath/utils/PolConverter.h
./askap/scimath/utils/PaddingUtils.tcc
./askap/scimath/utils/PhaseUnwrapper.tcc
./askap/scimath/utils/MultiDimArrayPlaneIter.tcc
./askap/scimath/utils/CasaBlobUtils.tcc
./askap/scimath/utils/FixedSizeCache.tcc

DESTINATION include/askap/scimath/utils
)

install (FILES

./askap/scimath/fft/FFTWrapper.h
DESTINATION include/askap/scimath/fft 

)
install (FILES
	
./askap/scimath/lsqr_solver/GlobalTypedefs.h
./askap/scimath/lsqr_solver/LSQRSolver.h
./askap/scimath/lsqr_solver/MathUtils.h
./askap/scimath/lsqr_solver/ModelDamping.h
./askap/scimath/lsqr_solver/ParallelTools.h
./askap/scimath/lsqr_solver/SparseMatrix.h

DESTINATION include/askap/scimath/lqsr_solver
)

target_link_libraries(
	askap_scimath
	askap::askap
	askap::parallel
	askap::imagemath
	lofar::Common
	lofar::Blob
	${CASACORE_LIBRARIES}
	${log4cxx_LIBRARY}
	${Boost_LIBRARIES}
	${FFTW_LIBRARIES}
	${GSL_LIBRARIES}
)
if (MPI_FOUND)
	target_link_libraries(askap_scimath
		${MPI_LIBRARIES}
	)
endif (MPI_FOUND)

target_include_directories(askap_scimath PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/askap/scimath>
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/askap/scimath>
)

install (
TARGETS askap_scimath
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

if (CPPUNIT_FOUND)
  
  include(CTest)
  enable_testing()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests/fft)
  add_subdirectory(tests/fft)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests/fitting)
  add_subdirectory(tests/fitting)
  
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests/utils)
  add_subdirectory(tests/utils)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests/lsqr_solver)
  add_subdirectory(tests/lsqr_solver)

endif()

