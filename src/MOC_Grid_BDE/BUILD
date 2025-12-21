cc_library(
    name = "MOC_Grid_BDE",

    srcs=[
        "MOC_GridCalc_BDE.cpp",
        "MOC_GridCalc_BDE_IO.cpp",
    ],
    hdrs=[
        "engineering_constants.hpp",
        "MOC_GridCalc_BDE.h",
        "dummyStruct.h",
        "stubs.h",
    ],
    copts = ["-Wno-everything"],
    #defines = ["LEGACY_CODE"],
    visibility = ["//visibility:public"],
    includes= ["src/MOC_Grid_BDE"],
    include_prefix = "MOC_Grid_BDE",
    strip_include_prefix = "",
)
