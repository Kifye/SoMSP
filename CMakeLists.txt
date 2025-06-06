cmake_minimum_required(VERSION 3.6)
project(Runnable VERSION 0.1.0)
enable_language (Fortran)

# set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran-11")
set (CMAKE_Fortran_FLAGS_RELEASE "-Ofast -flto -fmax-stack-var-size=65536")
set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check -fmax-stack-var-size=65536")
# set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fsanitize=address")

set(CMAKE_Fortran_MODULE_DIRECTORY "./modfiles")

add_library(regime src/regime.F90)

set(SOURCES
    src/constants.F90
    src/logging.f90
    src/functions/complex_oblate_swf.f90
    src/functions/complex_prolate_swf.f90
    src/functions/spheroidal.f90
    src/functions/legendre_functions.f90
    src/functions/param.f90
    src/integrals.f90
    src/mode_functors.f90
    src/matrix.f90
    src/spheroidal_tmatrix.f90
    src/tmatrix_conversion.f90
    src/scatterer.f90
    src/spheroidal_indicatrix.f90
    src/angle.f90
    src/contexts/contexts.f90
    src/wavelength_point.f90
    src/scattering_calculation.f90
    src/utils.f90
    src/spheroidal_scattering.f90
    src/spheroidal_initial.f90
    src/functions/legendre_functions.f90
    src/directions.f90
    src/calculation_model.f90
    src/tmatrix_extraction.f90
    src/geometry.f90
    src/contexts/spheroidal_context.f90
    src/contexts/scattering_context.f90
    src/contexts/far_context.f90
    src/contexts/barber_context.f90
    src/functions/bessel_function.f90
)
add_executable(main src/factors.f90 ${SOURCES})
add_executable(compare_spvv8 src/testing/compare_spvv8.f90 ${SOURCES})
add_executable(compare_spvv8_many_lambda src/testing/compare_spvv8_many_lambda.f90 ${SOURCES})
add_executable(check_pis src/testing/check_pis.f90 ${SOURCES})

if(USE_DOUBLE_PRECISION)
    message("utilize double precision")
    add_compile_definitions(USE_DOUBLE_PRECISION=True) 
    # target_compile_definitions(main PUBLIC -DUSE_DOUBLE_PRECISION=True) 
endif(USE_DOUBLE_PRECISION)

if(DEFINED LOG)
    message("print log to scattering.log")
    target_compile_definitions(main PUBLIC -DNEED_LOG) 
    if ("tmatrix" IN_LIST LOG)
        message("print tmatrices")
        target_compile_definitions(main PUBLIC -DNEED_LOG_TMATRIX) 
    endif()
    if ("time" IN_LIST LOG)
        message("print times")
        target_compile_definitions(main PUBLIC -DNEED_LOG_TIME) 
    endif()
    if ("amplitude" IN_LIST LOG)
        message("print amplitude")
        target_compile_definitions(main PUBLIC -DNEED_LOG_AMPLITUDE) 
    endif()
endif()
target_link_libraries(main regime)

add_library(spheroidal_scattering ${SOURCES})
target_link_libraries(spheroidal_scattering regime)

