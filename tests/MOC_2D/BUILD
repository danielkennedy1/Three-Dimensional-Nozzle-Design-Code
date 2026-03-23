cc_test(
    name = "moc_2d_test",
    srcs = ["moc_2d_test.cpp"],
    deps = [
        "//tests/fixtures:solver",
        "//src/MOC_2D/solver:solver",
        "//src/MOC_2D/common:common",
        "//src/MOC_2D/config:config",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
    copts = ["-std=c++23"],
    size = "medium",
)
