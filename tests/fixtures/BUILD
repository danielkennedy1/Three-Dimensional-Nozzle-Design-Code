cc_library(
    name = "legacy_2d",
    hdrs = ["legacy_2d.h"],
    srcs = ["legacy_2d.cpp"],
    deps = [
        "//src/MOC_Grid_BDE:MOC_Grid_BDE",
        "@googletest//:gtest",
    ],
    visibility = ["//tests:__subpackages__"],
)

cc_library(
    name = "solver",
    hdrs = ["solver.h"],
    srcs = ["solver.cpp"],
    deps = [
        "//src/MOC_2D/solver:solver",
        "@googletest//:gtest",
    ],
    visibility = ["//tests:__subpackages__"],
)
