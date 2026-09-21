add_library(mph STATIC
	mphcore_mod.F90
	mph_mod.F90
	mph_utils_mod.F90
	mph_props_mod.F90
	mph_test_mod.F90
	mph_plic_mod.F90
	mph_vof_mod.F90
	mph_pois_mod.F90
)

target_include_directories(mph
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(mph
    PRIVATE core
    PRIVATE ib
    PRIVATE flow
    PRIVATE ${MPI_Fortran_LIBRARIES}
)
