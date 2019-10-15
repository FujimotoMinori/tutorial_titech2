# Generated by lcg_generate_env...
if [ -z "${LCG_RELEASE_BASE}" ]; then
   export LCG_RELEASE_BASE=/cvmfs/atlas.cern.ch/repo/sw/software/21.2/sw/lcg/releases
fi
if [ -z "${LCG_PLATFORM}" ]; then
   export LCG_PLATFORM=x86_64-slc6-gcc62-opt
fi
export ROOTSYS=${LCG_RELEASE_BASE}/LCG_94a/ROOT/6.14.08/x86_64-slc6-gcc62-opt
export PYTHONHOME=${LCG_RELEASE_BASE}/LCG_94a/Python/2.7.15/x86_64-slc6-gcc62-opt
if [ -z "${GAUDI_ROOT}" ]; then
   export GAUDI_ROOT=/cvmfs/atlas.cern.ch/repo/sw/software/21.2/GAUDI/21.2.74/InstallArea/x86_64-slc6-gcc62-opt
fi
if [ -z "${JOBOPTSEARCHPATH}" ]; then
   export JOBOPTSEARCHPATH=${GAUDI_ROOT}/jobOptions
else
   export JOBOPTSEARCHPATH=${GAUDI_ROOT}/jobOptions:${JOBOPTSEARCHPATH}
fi
if [ -z "${DATAPATH}" ]; then
   export DATAPATH=${GAUDI_ROOT}/share
else
   export DATAPATH=${GAUDI_ROOT}/share:${DATAPATH}
fi
if [ -z "${PYTHONPATH}" ]; then
   export PYTHONPATH=${ROOTSYS}/lib:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysisExternals/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/python:${GAUDI_ROOT}/python:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysis/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/python
else
   export PYTHONPATH=${ROOTSYS}/lib:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysisExternals/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/python:${GAUDI_ROOT}/python:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysis/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/python:${PYTHONPATH}
fi
if [ -z "${PATH}" ]; then
   export PATH=${ROOTSYS}/bin:${PYTHONHOME}/bin:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysisExternals/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/bin:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysis/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/bin:${LCG_RELEASE_BASE}/LCG_94a/cairo/1.15.8/x86_64-slc6-gcc62-opt/bin:${LCG_RELEASE_BASE}/LCG_94a/graphviz/2.28.0/x86_64-slc6-gcc62-opt/bin:${LCG_RELEASE_BASE}/LCG_94a/doxygen/1.8.11/x86_64-slc6-gcc62-opt/bin
else
   export PATH=${ROOTSYS}/bin:${PYTHONHOME}/bin:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysisExternals/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/bin:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysis/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/bin:${LCG_RELEASE_BASE}/LCG_94a/cairo/1.15.8/x86_64-slc6-gcc62-opt/bin:${LCG_RELEASE_BASE}/LCG_94a/graphviz/2.28.0/x86_64-slc6-gcc62-opt/bin:${LCG_RELEASE_BASE}/LCG_94a/doxygen/1.8.11/x86_64-slc6-gcc62-opt/bin:${PATH}
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then
   export LD_LIBRARY_PATH=${LCG_RELEASE_BASE}/LCG_94a/Boost/1.66.0/x86_64-slc6-gcc62-opt/lib::${LCG_RELEASE_BASE}/LCG_94a/Boost/1.66.0/x86_64-slc6-gcc62-opt/lib::::${LCG_RELEASE_BASE}/LCG_94a/tbb/2018_U1/x86_64-slc6-gcc62-opt/lib:${PYTHONHOME}/lib:${ROOTSYS}/lib:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysisExternals/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/lib:${GAUDI_ROOT}/lib:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysis/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/lib:${LCG_RELEASE_BASE}/LCG_94a/cairo/1.15.8/x86_64-slc6-gcc62-opt/lib:${LCG_RELEASE_BASE}/LCG_94a/pango/1.40.13/x86_64-slc6-gcc62-opt/lib:${LCG_RELEASE_BASE}/LCG_94a/graphviz/2.28.0/x86_64-slc6-gcc62-opt/lib
else
   export LD_LIBRARY_PATH=${LCG_RELEASE_BASE}/LCG_94a/Boost/1.66.0/x86_64-slc6-gcc62-opt/lib::${LCG_RELEASE_BASE}/LCG_94a/Boost/1.66.0/x86_64-slc6-gcc62-opt/lib::::${LCG_RELEASE_BASE}/LCG_94a/tbb/2018_U1/x86_64-slc6-gcc62-opt/lib:${PYTHONHOME}/lib:${ROOTSYS}/lib:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysisExternals/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/lib:${GAUDI_ROOT}/lib:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysis/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/lib:${LCG_RELEASE_BASE}/LCG_94a/cairo/1.15.8/x86_64-slc6-gcc62-opt/lib:${LCG_RELEASE_BASE}/LCG_94a/pango/1.40.13/x86_64-slc6-gcc62-opt/lib:${LCG_RELEASE_BASE}/LCG_94a/graphviz/2.28.0/x86_64-slc6-gcc62-opt/lib:${LD_LIBRARY_PATH}
fi
if [ -z "${ROOT_INCLUDE_PATH}" ]; then
   export ROOT_INCLUDE_PATH=${LCG_RELEASE_BASE}/LCG_94a/Boost/1.66.0/x86_64-slc6-gcc62-opt/include::${LCG_RELEASE_BASE}/LCG_94a/Boost/1.66.0/x86_64-slc6-gcc62-opt/include::::${LCG_RELEASE_BASE}/LCG_94a/tbb/2018_U1/x86_64-slc6-gcc62-opt/include::${PYTHONHOME}/include/python2.7::${ROOTSYS}/include::/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysisExternals/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/include:${GAUDI_ROOT}/include:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysis/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/include:${LCG_RELEASE_BASE}/LCG_94a/cairo/1.15.8/x86_64-slc6-gcc62-opt/include::${LCG_RELEASE_BASE}/LCG_94a/pango/1.40.13/x86_64-slc6-gcc62-opt/include/pango-1.0::${LCG_RELEASE_BASE}/LCG_94a/graphviz/2.28.0/x86_64-slc6-gcc62-opt/include:
else
   export ROOT_INCLUDE_PATH=${LCG_RELEASE_BASE}/LCG_94a/Boost/1.66.0/x86_64-slc6-gcc62-opt/include::${LCG_RELEASE_BASE}/LCG_94a/Boost/1.66.0/x86_64-slc6-gcc62-opt/include::::${LCG_RELEASE_BASE}/LCG_94a/tbb/2018_U1/x86_64-slc6-gcc62-opt/include::${PYTHONHOME}/include/python2.7::${ROOTSYS}/include::/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysisExternals/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/include:${GAUDI_ROOT}/include:/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AthAnalysis/21.2.74/InstallArea/x86_64-slc6-gcc62-opt/include:${LCG_RELEASE_BASE}/LCG_94a/cairo/1.15.8/x86_64-slc6-gcc62-opt/include::${LCG_RELEASE_BASE}/LCG_94a/pango/1.40.13/x86_64-slc6-gcc62-opt/include/pango-1.0::${LCG_RELEASE_BASE}/LCG_94a/graphviz/2.28.0/x86_64-slc6-gcc62-opt/include::${ROOT_INCLUDE_PATH}
fi
