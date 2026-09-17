add_library(multiphase STATIC
	mphcore_mod.F90
	mph_mod.F90
	mph_utils_mod.F90
	mph_props_mod.F90
	mph_test_mod.F90
	mph_plic_mod.F90
	mph_vof_mod.F90
)

target_include_directories(multiphase
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(multiphase
    PRIVATE core
    PRIVATE ib
    PRIVATE flow
    PRIVATE ${MPI_Fortran_LIBRARIES}
)
