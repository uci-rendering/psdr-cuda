#
#  Copyright (c) 2018 NVIDIA Corporation.  All rights reserved.
#
#  NVIDIA Corporation and its licensors retain all intellectual property and proprietary
#  rights in and to this software, related documentation and any modifications thereto.
#  Any use, reproduction, disclosure or distribution of this software and related
#  documentation without an express license agreement from NVIDIA Corporation is strictly
#  prohibited.
#
#  TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, THIS SOFTWARE IS PROVIDED *AS IS*
#  AND NVIDIA AND ITS SUPPLIERS DISCLAIM ALL WARRANTIES, EITHER EXPRESS OR IMPLIED,
#  INCLUDING, BUT NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#  PARTICULAR PURPOSE.  IN NO EVENT SHALL NVIDIA OR ITS SUPPLIERS BE LIABLE FOR ANY
#  SPECIAL, INCIDENTAL, INDIRECT, OR CONSEQUENTIAL DAMAGES WHATSOEVER (INCLUDING, WITHOUT
#  LIMITATION, DAMAGES FOR LOSS OF BUSINESS PROFITS, BUSINESS INTERRUPTION, LOSS OF
#  BUSINESS INFORMATION, OR ANY OTHER PECUNIARY LOSS) ARISING OUT OF THE USE OF OR
#  INABILITY TO USE THIS SOFTWARE, EVEN IF NVIDIA HAS BEEN ADVISED OF THE POSSIBILITY OF
#  SUCH DAMAGES
#

# Find OpenEXR package.

# Optional input variable: OpenEXR_ROOT
# Output variables:
#   OpenEXR_FOUND
#   OpenEXR_INCLUDE_DIR
#   OpenEXR_LIBRARIES
#   OpenEXR_VERSION

set( OpenEXR_LIB_NAMES IlmImf Half Iex Imath IlmThread )

# If OpenEXR_ROOT has changed, unset variables that depend upon it.
set( OpenEXR_ROOT "" CACHE PATH "Path to OpenEXR installation directory" )
if( NOT "${OpenEXR_ROOT}" STREQUAL "${OpenEXR_ROOT_PREVIOUS}" )
  message( "New value detected for OpenEXR_ROOT: ${OpenEXR_ROOT}" )
  unset( OpenEXR_INCLUDE_DIR )
  unset( OpenEXR_LIBRARIES )
  unset( OpenEXR_LIB_DIR )
  foreach( LIB ${OpenEXR_LIB_NAMES} )
    unset( OpenEXR_${LIB}_RELEASE )
    unset( OpenEXR_${LIB}_DEBUG )
  endforeach( LIB )
  unset( OpenEXR_VERSION )
endif()
set( OpenEXR_ROOT_PREVIOUS "${OpenEXR_ROOT}" CACHE PATH "Previous path to OpenEXR" FORCE )

# Find OpenEXR includes.
find_path( OpenEXR_INCLUDE_DIR ImfOutputFile.h
	   HINTS "${OpenEXR_ROOT}/include/OpenEXR" )
mark_as_advanced( OpenEXR_INCLUDE_DIR )

# Get version number from header, which we need for the library names.
set( OpenEXR_VERSION "" CACHE STRING "OpenEXR version string" )
set( CONFIG_H "${OpenEXR_INCLUDE_DIR}/OpenEXRConfig.h" )
if( NOT OpenEXR_VERSION AND EXISTS "${CONFIG_H}" )
  message( "Reading OpenEXR version from ${CONFIG_H}" )
  file( STRINGS "${CONFIG_H}" VERSION_STRING
    REGEX "#define OPENEXR_VERSION_STRING" )
  string( REGEX REPLACE ".*\"([0-9.]+)\".*" "\\1" VERSION_STRING "${VERSION_STRING}" )
  set( OpenEXR_VERSION "${VERSION_STRING}" CACHE STRING "OpenEXR version string" FORCE )
endif()
string( REGEX REPLACE "^([0-9]+).*" "\\1" VERSION_MAJOR "${OpenEXR_VERSION}" )
string( REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${OpenEXR_VERSION}" )
set( VERSION_SUFFIX "${VERSION_MAJOR}_${VERSION_MINOR}" )

# Allow location of library directory to be overridden.
set( OpenEXR_LIB_DIR "${OpenEXR_ROOT}/lib" CACHE PATH "Path to OpenEXR libraries" )
mark_as_advanced( OpenEXR_LIB_DIR )

# Find OpenEXR libraries.
set( OpenEXR_LIBRARIES "" )
foreach( LIB ${OpenEXR_LIB_NAMES} )
  find_library( OpenEXR_${LIB}_RELEASE
  		NAMES "${LIB}_s" "${LIB}-${VERSION_SUFFIX}_s" "${LIB}"
 		HINTS "${OpenEXR_LIB_DIR}" )
  mark_as_advanced( OpenEXR_${LIB}_RELEASE )
  if( OpenEXR_${LIB}_RELEASE )
    list( APPEND OpenEXR_LIBRARIES optimized "${OpenEXR_${LIB}_RELEASE}" )
  endif()

  find_library( OpenEXR_${LIB}_DEBUG
  		NAMES "${LIB}_s_d" "${LIB}-${VERSION_SUFFIX}_s_d"
 		HINTS "${OpenEXR_LIB_DIR}" )
  mark_as_advanced( OpenEXR_${LIB}_DEBUG )
  if( OpenEXR_${LIB}_DEBUG )
    list( APPEND OpenEXR_LIBRARIES debug "${OpenEXR_${LIB}_DEBUG}" )
  elseif( OpenEXR_${LIB}_RELEASE )
    # Fallback: use release libraries if no debug libraries were found.
    list( APPEND OpenEXR_LIBRARIES debug "${OpenEXR_${LIB}_RELEASE}" )
  endif()
endforeach( LIB )

include( FindPackageHandleStandardArgs )

# find_package_handle_standard_args reports the value of the first variable
# on success, so make sure this is the actual OpenEXR library
find_package_handle_standard_args( OpenEXR
  REQUIRED_VARS
    OpenEXR_IlmImf_RELEASE OpenEXR_Half_RELEASE OpenEXR_Iex_RELEASE OpenEXR_Imath_RELEASE OpenEXR_IlmThread_RELEASE
    OpenEXR_INCLUDE_DIR
  VERSION_VAR OpenEXR_VERSION )
