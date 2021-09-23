# This file is part of https://github.com/PalamaraLab/ASMC which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# If a VCPKG toolchain file is not defined, but the expected file exists, use it
if (NOT EXISTS ${CMAKE_SOURCE_DIR}/DataModule/README.md)
    message(FATAL_ERROR "
The data module ${ASMC_data_module_dir} does not exist, and it is required for ASMC.
Please either get all submodules when you clone ASMC:
$ git clone --recurse-submodules https://github.com/PalamaraLab/ASMC.git
or, at minimum, initialise the data module. From the ASMC directory:
$ git submodule update --init DataModule
Please see this quickstart guide for further information:
https://github.com/PalamaraLab/ASMC/blob/main/docs/quickstart_user.md
")
endif ()
