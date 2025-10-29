add_executable(
    wline
    main.cpp
    wilsonline.cpp
    ipglasma.cpp
)

target_link_libraries(
    wline
    PRIVATE
    GSL::gsl
    GSL::gslcblas
)

set_target_properties(wline PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

