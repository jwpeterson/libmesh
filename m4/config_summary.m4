# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   git log -n1 m4/config_summary.m4
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

######################################################################################
echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version.................... : $PACKAGE-$VERSION
echo
echo C++ compiler type.................. : $GXX_VERSION
echo C++ compiler....................... : $CXX
echo C compiler......................... : $CC
echo Fortran compiler................... : $FC
echo Build Methods...................... : $METHODS
echo " "
for method in ${METHODS}; do
     case "${method}" in
         opt)
echo CPPFLAGS...\(opt\)................... : $CPPFLAGS_OPT $CPFLAGS
echo CXXFLAGS...\(opt\)................... : $CXXFLAGS_OPT $CXXFLAGS
echo CFLAGS.....\(opt\)................... : $CFLAGS_OPT $CFLAGS
     ;;
         devel)
echo CPPFLAGS...\(devel\)................. : $CPPFLAGS_DEVEL $CPFLAGS
echo CXXFLAGS...\(devel\)................. : $CXXFLAGS_DEVEL $CXXFLAGS
echo CFLAGS.....\(devel\)................. : $CFLAGS_DEVEL $CFLAGS
     ;;
         dbg)
echo CPPFLAGS...\(dbg\)................... : $CPPFLAGS_DBG $CPFLAGS
echo CXXFLAGS...\(dbg\)................... : $CXXFLAGS_DBG $CXXFLAGS
echo CFLAGS.....\(dbg\)................... : $CFLAGS_DBG $CFLAGS
     ;;
         prof)
echo CPPFLAGS...\(prof\).................. : $CPPFLAGS_PROF $CPFLAGS
echo CXXFLAGS...\(prof\).................. : $CXXFLAGS_PROF $CXXFLAGS
echo CFLAGS.....\(prof\).................. : $CFLAGS_PROF $CFLAGS
     ;;
         oprof)
echo CPPFLAGS...\(oprof\)................. : $CPPFLAGS_OPROF $CPFLAGS
echo CXXFLAGS...\(oprof\)................. : $CXXFLAGS_OPROF $CXXFLAGS
echo CFLAGS.....\(oprof\)................. : $CFLAGS_OPROF $CFLAGS
     esac
     echo " "
done
echo Install dir........................ : $prefix
echo Build user......................... : $USER
echo Build host......................... : $BUILD_HOST
echo Build architecture................. : $BUILD_ARCH
echo Git revision....................... : $BUILD_VERSION

######################################################################################
echo
echo Library Features:
echo '  library warnings................. :' $enablewarnings
echo '  library deprecated code support.. :' $enabledeprecated
echo '  adaptive mesh refinement......... :' $enableamr
echo '  blocked matrix/vector storage.... :' $enableblockedstorage
echo '  complex variables................ :' $enablecomplex
echo '  example suite.................... :' $enableexamples
echo '  ghosted vectors.................. :' $enableghosted
echo '  high-order shape functions....... :' $enablepfem
echo '  unique-id support................ :' $enableuniqueid
echo '  id size (boundaries)............. :' $boundary_bytes bytes
echo '  id size (dofs)................... :' $dof_bytes bytes
AS_IF([test "x$enableuniqueid" = "xyes"],
      [echo '  id size (unique)................. :' $unique_bytes bytes])
echo '  id size (processors)............. :' $processor_bytes bytes
echo '  id size (subdomains)............. :' $subdomain_bytes bytes
echo '  infinite elements................ :' $enableifem
echo '  Dirichlet constraints............ :' $enabledirichlet
echo '  node constraints................. :' $enablenodeconstraint
echo '  parallel mesh.................... :' $enableparmesh
echo '  performance logging.............. :' $enableperflog
echo '  periodic boundary conditions..... :' $enableperiodic
echo '  reference counting............... :' $enablerefct
echo '  shape function 2nd derivatives... :' $enablesecond
echo '  stack trace files................ :' $enabletracefiles
echo '  track node valence............... :' $enablenodevalence
echo '  variational smoother............. :' $enablevsmoother
echo '  xdr binary I/O................... :' $enablexdr
if (test "x$enablelegacyincludepaths" = "xyes"); then
echo '  non-prefixed include paths....... :' $enablelegacyincludepaths ***LEGACY FEATURE***
fi
if (test "x$enablelegacyusingnamespace" = "xyes"); then
echo '  adding using namespace libMesh... :' $enablelegacyusingnamespace ***LEGACY FEATURE***
fi
if (test "x$enabledefaultcommworld" = "xyes"); then
echo '  providing libMesh::CommWorld..... :' $enabledefaultcommworld ***LEGACY FEATURE***
fi



######################################################################################
if (test "x$enableoptional" = "xyes"); then
  echo
  echo Optional Packages:
  echo '  'boost............................ : $enableboost
  echo '  'capnproto........................ : $enablecapnproto
  echo '  'cppunit.......................... : $enablecppunit
  echo '  'curl............................. : $enablecurl
  echo '  'eigen............................ : $enableeigen
  echo '  'exodus........................... : $enableexodus
  AS_IF([test "x$exodusversion" != "xno"],
        [echo '     'version....................... : $exodusversion])
  echo '  'fparser.......................... : $enablefparser
  AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdevel" = "xno"],
        [echo '     'build from version............ : release])
  AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdevel" = "xyes"],
        [echo '     'build from version............ : devel])
  AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdebugging" = "xyes"],
        [echo '     'fparser debugging............. : enabled])
  echo '  'glpk............................. : $enableglpk
  echo '  'gmv.............................. : $enablegmv
  echo '  'gzstream......................... : $enablegz
  echo '  'hdf5............................. : $enablehdf5
  echo '  'laspack.......................... : $enablelaspack
  echo '  'libhilbert....................... : $enablelibhilbert
  echo '  'metis............................ : $enablemetis
  echo '  'mpi.............................. : $enablempi
  echo '  'nanoflann........................ : $enablenanoflann
  echo '  'nemesis.......................... : $enablenemesis
  AS_IF([test "x$nemesisversion" != "xno"],
        [echo '     'version....................... : $nemesisversion])
  echo '  'netcdf........................... : $enablenetcdf
  AS_IF([test "x$netcdfversion" != "xno"],
        [echo '     'version....................... : $netcdfversion])
  echo '  'nlopt............................ : $enablenlopt
  echo '  'parmetis......................... : $enableparmetis
  echo '  'petsc............................ : $enablepetsc
  AS_IF([test "x$enablepetsc" = "xyes"],
        [echo '     'version....................... : $petscversion])
  echo '  'qhull............................ : $enableqhull
  echo '  'sfcurves......................... : $enablesfc
  echo '  'slepc............................ : $enableslepc
  AS_IF([test "x$enableslepc" = "xyes"],
        [echo '     'version....................... : $slepcversion])
  echo '  'thread model..................... : $found_thread_model
  echo '  'c++ rtti ........................ : $ac_cv_cxx_rtti
  echo '  'tecio............................ : $enabletecio
  echo '  'tecplot...\(vendor binaries\)...... : $enabletecplot
  echo '  'tetgen........................... : $enabletetgen
  echo '  'triangle......................... : $enabletriangle
  echo '  'trilinos......................... : $enabletrilinos
  AS_IF([test "x$enabletrilinos" = "xyes"],
        [
          echo '     'AztecOO....................... : $enableaztecoo
          echo '     'NOX........................... : $enablenox
          echo '     'ML............................ : $enableml
          echo '     'Tpetra........................ : $enabletpetra
          echo '     'DTK........................... : $enabledtk
          echo '     'Ifpack........................ : $enableifpack
          echo '     'Epetra........................ : $enableepetra
          echo '     'EpetraExt..................... : $enableepetraext
        ])
  echo '  'vtk.............................. : $enablevtk
  AS_IF([test "x$enablevtk" = "xyes"],
        [echo '     'version....................... : $vtkversion])
  dnl blank line
  echo
  AS_IF([test "x$libmesh_optional_INCLUDES" != "x"],
        [
          echo '  'libmesh_optional_INCLUDES........ : $libmesh_optional_INCLUDES
          echo
        ])
  AS_IF([test "x$libmesh_optional_LIBS" != "x"],
        [
          echo '  'libmesh_optional_LIBS............ : $libmesh_optional_LIBS
          echo
        ])
fi
echo '-------------------------------------------------------------------------------'

echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
