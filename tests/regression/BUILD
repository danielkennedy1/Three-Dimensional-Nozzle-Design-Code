cc_test(
    name = "moc_regression_test",
    srcs = ["moc_regression_test.cpp"],
    deps = [
        "//tests/fixtures:legacy_2d",
        "//tests/fixtures:solver",
        "//src/MOC_2D/solver:solver",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
    copts = ["-std=c++23"],
    size = "medium",
)
