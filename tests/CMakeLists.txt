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


# 2. Unit tests: Response_Tabulation
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

# 3. Unit tests: Special_Functions
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

# 4. Unit tests: Wavefunctions
add_executable(test_Wavefunctions
	test_Wavefunctions.cpp)

target_link_libraries(test_Wavefunctions 
	PRIVATE
		libdarkarc
		gtest_main ) #contains the main function

target_include_directories(test_Wavefunctions
	PRIVATE
		${GENERATED_DIR} )

target_compile_options(test_Wavefunctions PUBLIC -Wall -pedantic)

install(TARGETS test_Wavefunctions DESTINATION ${TESTS_DIR})

add_test(NAME test_Wavefunctions COMMAND test_Wavefunctions
	WORKING_DIRECTORY ${TESTS_DIR} )