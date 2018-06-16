add_executable(
    tests
    main.cpp
    distribution_test.cpp
    bit_reversal_test.cpp
    )

target_link_libraries(
    tests
    util
    )
