# CMakeLists.txt for the tests directory

# 1. Define the test executable and its source files
add_executable(runUnitTests
    test_Logger.cpp
)

# 2. Link the test executable against our toolbox library
target_link_libraries(runUnitTests
    PRIVATE
        IsoToolbox
)

# 3. Register the executable as a test with CTest.
#    NAME sets the name of the test, COMMAND is the program to run.
add_test(
    NAME LoggerTest
    COMMAND runUnitTests
)