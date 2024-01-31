# Install script for directory: /home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local/HDF_Group/HDF5/1.10.5")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/home/rfaure/miniconda3/bin/x86_64-conda-linux-gnu-objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/hdf5.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5api_adpt.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5public.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Apublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5ACpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Cpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Dpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Epubgen.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Epublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Fpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDcore.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDdirect.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDfamily.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDlog.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDmpi.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDmpio.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDmulti.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDsec2.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDstdio.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5FDwindows.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Gpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Ipublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Lpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5MMpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Opublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Ppublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5PLextern.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5PLpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Rpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Spublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Tpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Zpublic.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5Epubgen.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5version.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5overflow.h;/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5/H5pubconf.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/include/hdf5" TYPE FILE FILES
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/hdf5.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5api_adpt.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5public.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Apublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5ACpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Cpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Dpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Epubgen.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Epublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Fpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDcore.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDdirect.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDfamily.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDlog.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDmpi.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDmpio.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDmulti.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDsec2.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDstdio.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDwindows.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Gpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Ipublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Lpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5MMpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Opublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Ppublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5PLextern.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5PLpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Rpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Spublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Tpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Zpublic.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5Epubgen.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5version.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/gatb-core/gatb-core/thirdparty/hdf5/src/H5overflow.h"
    "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty/hdf5/H5pubconf.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty/hdf5/CMakeFiles/hdf5-1.10.5.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/home/rfaure/Documents/these/MSR/wonderasm/bcalm/build/ext/gatb-core/thirdparty/hdf5/CMakeFiles/h5cc")
endif()

