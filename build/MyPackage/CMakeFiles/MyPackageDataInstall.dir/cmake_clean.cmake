file(REMOVE_RECURSE
  "../x86_64-centos7-gcc62-opt/include/MyPackage"
  "../x86_64-centos7-gcc62-opt/python/MyPackage/__init__.py"
  "../x86_64-centos7-gcc62-opt/jobOptions/MyPackage/MyPackageAlgJobOptions.py"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/MyPackageDataInstall.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
