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
	file(GLOB FILES "*.cpp")
	list(REMOVE_ITEM FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
	
	add_library(libdarkart STATIC
    ${FILES})

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