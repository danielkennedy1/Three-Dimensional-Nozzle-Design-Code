cc_library(
    name = "MOC_2D",
    srcs=[
        "MOC_GridCalc.cpp",
        "MOC_GridCalc_IO.cpp",
    ],
    hdrs=[
        "engineering_constants.hpp",
        "MOC_GridCalc.h",
        "dummyStruct.h",
        "stubs.h",
    ],
    copts = ["-Wno-everything"],
    visibility = ["//tests:__subpackages__"],
    includes= ["src/MOC_2D"],
    include_prefix = "MOC_2D",
    strip_include_prefix = "",
)
