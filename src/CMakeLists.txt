# OPTION A: Build an executable and a static library
	# Executable
	add_executable(DarkARC++
		main.cpp)

	target_compile_options(DarkARC++ PUBLIC -Wall -pedantic)

	target_link_libraries(DarkARC++
		PUBLIC
			coverage_config
			libdarkarc)

	target_include_directories(DarkARC++
		PRIVATE
			${GENERATED_DIR} )

	install(TARGETS DarkARC++ DESTINATION ${BIN_DIR})

	# Static library
	add_library(libdarkarc STATIC
	    Atomic_Responses.cpp
	    Configuration.cpp
	    Response_Tabulation.cpp
	    Special_Functions.cpp
	    Wavefunctions.cpp)

	set_target_properties(libdarkarc PROPERTIES PREFIX "")

	target_compile_options(libdarkarc PUBLIC -Wall -pedantic)

	target_link_libraries(libdarkarc
		PRIVATE
			${Boost_LIBRARIES}
			GSL::gsl
			GSL::gslcblas
			"${ARB}"
			OpenMP::OpenMP_CXX
		PUBLIC
			coverage_config
			libphysica)

	target_include_directories(libdarkarc
		PRIVATE
			${GENERATED_DIR}
		PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}
			${INCLUDE_DIR})

	install(TARGETS libdarkarc DESTINATION ${LIB_DIR})