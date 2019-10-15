#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "WorkDir::MyPackageLib" for configuration "RelWithDebInfo"
set_property(TARGET WorkDir::MyPackageLib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(WorkDir::MyPackageLib PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libMyPackageLib.so"
  IMPORTED_SONAME_RELWITHDEBINFO "libMyPackageLib.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS WorkDir::MyPackageLib )
list(APPEND _IMPORT_CHECK_FILES_FOR_WorkDir::MyPackageLib "${_IMPORT_PREFIX}/lib/libMyPackageLib.so" )

# Import target "WorkDir::MyPackage" for configuration "RelWithDebInfo"
set_property(TARGET WorkDir::MyPackage APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(WorkDir::MyPackage PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libMyPackage.so"
  IMPORTED_NO_SONAME_RELWITHDEBINFO "TRUE"
  )

list(APPEND _IMPORT_CHECK_TARGETS WorkDir::MyPackage )
list(APPEND _IMPORT_CHECK_FILES_FOR_WorkDir::MyPackage "${_IMPORT_PREFIX}/lib/libMyPackage.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
