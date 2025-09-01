# Define the test executable
add_executable(run_tests
    test_Logger.cpp
)

# Link tests against the libraries they are testing
target_link_libraries(run_tests
    PRIVATE
        IsoToolbox
)

# Add the test to be run by 'ctest'
add_test(NAME TestLogger COMMAND run_tests)