# Define the executable target
add_executable(process
    src/main.cpp
)

# Link the executable against all the libraries it needs to run
target_link_libraries(process
    PRIVATE
        # The executable is the final product, so its dependencies are private
        HistManager
        Selection
        DataModel
        IsoToolbox
)