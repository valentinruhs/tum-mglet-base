add_library(multiphase STATIC
	multiphasecore_mod.F90
	multiphase_mod.F90
)

target_include_directories(multiphase
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(multiphase
    PRIVATE core
    PRIVATE ib
    PRIVATE ${MPI_Fortran_LIBRARIES}
)
