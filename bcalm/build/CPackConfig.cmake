# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_BUILD_SOURCE_DIRS "/home/rfaure/Documents/these/MSR/wonderasm/bcalm;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENTS_ALL "Unspecified;headers;libraries")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/usr/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "bcalm built using CMake")
set(CPACK_GENERATOR "TGZ")
set(CPACK_GFORGE_PROJECT_NAME "gatb-tools")
set(CPACK_INSTALL_CMAKE_PROJECTS "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build;bcalm;ALL;/")
set(CPACK_INSTALL_PREFIX "/usr/local/HDF_Group/HDF5/1.10.5")
set(CPACK_MODULE_PATH "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/cmake;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/cmake")
set(CPACK_NSIS_DISPLAY_NAME "bcalm ")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_PACKAGE_NAME "bcalm ")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_OUTPUT_CONFIG_FILE "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/usr/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "bcalm built using CMake")
set(CPACK_PACKAGE_FILE_NAME "bcalm-binaries-v2.2.3-Linux")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "bcalm ")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "bcalm ")
set(CPACK_PACKAGE_NAME "bcalm")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "Humanity")
set(CPACK_PACKAGE_VERSION "")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_RESOURCE_FILE_LICENSE "/usr/share/cmake/Templates/CPack.GenericLicense.txt")
set(CPACK_RESOURCE_FILE_README "/usr/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_RESOURCE_FILE_WELCOME "/usr/share/cmake/Templates/CPack.GenericWelcome.txt")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "^/home/rfaure/Documents/these/MSR/wonderasm/bcalm/.git/;^/home/rfaure/Documents/these/MSR/wonderasm/bcalm/.gitmodules;^/home/rfaure/Documents/these/MSR/wonderasm/bcalm/.gitignore;^/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/;^/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/.cproject;^/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/.git/;^/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/.project;^/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/.gitignore")
set(CPACK_SOURCE_INSTALLED_DIRECTORIES "/home/rfaure/Documents/these/MSR/wonderasm/bcalm;.;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core;thirdparty/gatb-core;;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/cmake;cmake")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/CPackSourceConfig.cmake")
set(CPACK_SYSTEM_NAME "Linux")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "Linux")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
