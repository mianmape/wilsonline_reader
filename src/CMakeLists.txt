
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
