# OPTION A: Build an executable and a static library
	# Executable
	add_executable(DarkART
		main.cpp)

	target_compile_options(DarkART PUBLIC -Wall -pedantic)

	target_link_libraries(DarkART
		PUBLIC
			coverage_config
			libdarkart)

	target_include_directories(DarkART
		PRIVATE
			${GENERATED_DIR} )

	install(TARGETS DarkART DESTINATION ${BIN_DIR})

	# Static library
	add_library(libdarkart STATIC
	    Atomic_Responses.cpp
	    Configuration.cpp
	    Radial_Integrator.cpp
	    Response_Tabulation.cpp
	    Special_Functions.cpp
	    Wavefunctions_Final.cpp
	    Wavefunctions_Initial.cpp)

	set_target_properties(libdarkart PROPERTIES PREFIX "")

	target_compile_options(libdarkart PUBLIC -Wall -pedantic)

	target_link_libraries(libdarkart
		PRIVATE
			${Boost_LIBRARIES}
			GSL::gsl
			GSL::gslcblas
			"${ARB}"
			OpenMP::OpenMP_CXX
		PUBLIC
			coverage_config
			libphysica)

	target_include_directories(libdarkart
		PRIVATE
			${GENERATED_DIR}
		PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}
			${INCLUDE_DIR})

	install(TARGETS libdarkart DESTINATION ${LIB_DIR})