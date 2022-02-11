# 0. Test executable
add_test(NAME test_DarkARC_Executable COMMAND DarkARC++ .test.cfg
	WORKING_DIRECTORY ${BIN_DIR})

# 1. Unit tests: Atomic_Responses
add_executable(test_Atomic_Responses
	test_Atomic_Responses.cpp)

target_link_libraries(test_Atomic_Responses 
	PRIVATE
		libdarkarc
		gtest_main ) #contains the main function

target_include_directories(test_Atomic_Responses
	PRIVATE
		${GENERATED_DIR} )

target_compile_options(test_Atomic_Responses PUBLIC -Wall -pedantic)

install(TARGETS test_Atomic_Responses DESTINATION ${TESTS_DIR})

add_test(NAME test_Atomic_Responses COMMAND test_Atomic_Responses
	WORKING_DIRECTORY ${TESTS_DIR} )

# 2. Unit tests: Configuration
add_executable(test_Configuration
test_Configuration.cpp)

target_link_libraries(test_Configuration 
PRIVATE
	libdarkarc
	gtest_main ) #contains the main function

target_include_directories(test_Configuration
PRIVATE
	${GENERATED_DIR} )

target_compile_options(test_Configuration PUBLIC -Wall -pedantic)

install(TARGETS test_Configuration DESTINATION ${TESTS_DIR})

add_test(NAME test_Configuration COMMAND test_Configuration
WORKING_DIRECTORY ${TESTS_DIR} )


# 3. Unit tests: Radial_Integrator
add_executable(test_Radial_Integrator
	test_Radial_Integrator.cpp)

target_link_libraries(test_Radial_Integrator 
	PRIVATE
		libdarkarc
		gtest_main ) #contains the main function

target_include_directories(test_Radial_Integrator
	PRIVATE
		${GENERATED_DIR} )

target_compile_options(test_Radial_Integrator PUBLIC -Wall -pedantic)

install(TARGETS test_Radial_Integrator DESTINATION ${TESTS_DIR})

add_test(NAME test_Radial_Integrator COMMAND test_Radial_Integrator
	WORKING_DIRECTORY ${TESTS_DIR} )

# 4. Unit tests: Response_Tabulation
add_executable(test_Response_Tabulation
	test_Response_Tabulation.cpp)

target_link_libraries(test_Response_Tabulation 
	PRIVATE
		libdarkarc
		gtest_main ) #contains the main function

target_include_directories(test_Response_Tabulation
	PRIVATE
		${GENERATED_DIR} )

target_compile_options(test_Response_Tabulation PUBLIC -Wall -pedantic)

install(TARGETS test_Response_Tabulation DESTINATION ${TESTS_DIR})

add_test(NAME test_Response_Tabulation COMMAND test_Response_Tabulation
	WORKING_DIRECTORY ${TESTS_DIR} )

# 5. Unit tests: Special_Functions
add_executable(test_Special_Functions
test_Special_Functions.cpp)

target_link_libraries(test_Special_Functions 
PRIVATE
	libdarkarc
	gtest_main ) #contains the main function

target_include_directories(test_Special_Functions
PRIVATE
	${GENERATED_DIR} )

target_compile_options(test_Special_Functions PUBLIC -Wall -pedantic)

install(TARGETS test_Special_Functions DESTINATION ${TESTS_DIR})

add_test(NAME test_Special_Functions COMMAND test_Special_Functions
WORKING_DIRECTORY ${TESTS_DIR} )

# 6. Unit tests: Wavefunctions_Final
add_executable(test_Wavefunctions_Final
	test_Wavefunctions_Final.cpp)

target_link_libraries(test_Wavefunctions_Final 
	PRIVATE
		libdarkarc
		gtest_main ) #contains the main function

target_include_directories(test_Wavefunctions_Final
	PRIVATE
		${GENERATED_DIR} )

target_compile_options(test_Wavefunctions_Final PUBLIC -Wall -pedantic)

install(TARGETS test_Wavefunctions_Final DESTINATION ${TESTS_DIR})

add_test(NAME test_Wavefunctions_Final COMMAND test_Wavefunctions_Final
	WORKING_DIRECTORY ${TESTS_DIR} )

# 7. Unit tests: Wavefunctions_Initial
add_executable(test_Wavefunctions_Initial
	test_Wavefunctions_Initial.cpp)

target_link_libraries(test_Wavefunctions_Initial 
	PRIVATE
		libdarkarc
		gtest_main ) #contains the main function

target_include_directories(test_Wavefunctions_Initial
	PRIVATE
		${GENERATED_DIR} )

target_compile_options(test_Wavefunctions_Initial PUBLIC -Wall -pedantic)

install(TARGETS test_Wavefunctions_Initial DESTINATION ${TESTS_DIR})

add_test(NAME test_Wavefunctions_Initial COMMAND test_Wavefunctions_Initial
	WORKING_DIRECTORY ${TESTS_DIR} )