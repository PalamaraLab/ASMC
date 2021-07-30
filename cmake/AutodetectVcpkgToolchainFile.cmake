# This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# If a VCPKG toolchain file is not defined, but the expected file exists, use it
if (NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    if (EXISTS ${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake)
        set(vcpkg_toolchain_file ${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake)
        message(STATUS "Detected vcpkg toolchain file at ${vcpkg_toolchain_file}")
        set(CMAKE_TOOLCHAIN_FILE ${vcpkg_toolchain_file})
    endif ()
endif ()
