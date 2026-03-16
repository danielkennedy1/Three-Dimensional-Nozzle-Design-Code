cc_test(
    name = "moc_regression_test",
    srcs = ["moc_regression_test.cpp"],
    deps = [
        "//src/MOC_Grid_BDE:MOC_Grid_BDE",
        "//src/MOC_2D:MOC_2D",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
    copts = ["-std=c++23"],
    size = "medium",
)
