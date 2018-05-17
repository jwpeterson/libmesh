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
AS_ECHO([])
echo '----------------------------------- SUMMARY -----------------------------------'
AS_ECHO([])
echo Package version.................... : $PACKAGE-$VERSION
AS_ECHO([])
echo C++ compiler type.................. : $GXX_VERSION
echo C++ compiler....................... : $CXX
echo C compiler......................... : $CC
echo Fortran compiler................... : $FC
echo Build Methods...................... : $METHODS
AS_ECHO([])
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
     AS_ECHO([])
done
echo Install dir........................ : $prefix
echo Build user......................... : $USER
echo Build host......................... : $BUILD_HOST
echo Build architecture................. : $BUILD_ARCH
echo Git revision....................... : $BUILD_VERSION

######################################################################################
AS_ECHO([])
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
  AS_ECHO([])
  echo Optional Packages:
  AS_ECHO(["  boost............................ : $enableboost"])
  AS_ECHO(["  capnproto........................ : $enablecapnproto"])
  AS_ECHO(["  cppunit.......................... : $enablecppunit"])
  AS_ECHO(["  curl............................. : $enablecurl"])
  AS_ECHO(["  eigen............................ : $enableeigen"])
  AS_ECHO(["  exodus........................... : $enableexodus"])
  AS_IF([test "x$exodusversion" != "xno"],
        [echo '     'version....................... : $exodusversion])
  AS_ECHO(["  fparser.......................... : $enablefparser"])
  AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdevel" = "xno"],
        [AS_ECHO(["     build from version............ : release"])])
  AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdevel" = "xyes"],
        [AS_ECHO(["     build from version............ : devel"])])
  AS_IF([test "x$enablefparser" = "xyes" && test "x$enablefparserdebugging" = "xyes"],
        [AS_ECHO(["     fparser debugging............. : enabled"])])
  AS_ECHO(["  glpk............................. : $enableglpk"])
  AS_ECHO(["  gmv.............................. : $enablegmv"])
  AS_ECHO(["  gzstream......................... : $enablegz"])
  AS_ECHO(["  hdf5............................. : $enablehdf5"])
  AS_ECHO(["  laspack.......................... : $enablelaspack"])
  AS_ECHO(["  libhilbert....................... : $enablelibhilbert"])
  AS_ECHO(["  metis............................ : $enablemetis"])
  AS_ECHO(["  mpi.............................. : $enablempi"])
  AS_ECHO(["  nanoflann........................ : $enablenanoflann"])
  AS_ECHO(["  nemesis.......................... : $enablenemesis"])
  AS_IF([test "x$nemesisversion" != "xno"],
        [echo '     'version....................... : $nemesisversion])
  AS_ECHO(["  netcdf........................... : $enablenetcdf"])
  AS_IF([test "x$netcdfversion" != "xno"],
        [echo '     'version....................... : $netcdfversion])
  AS_ECHO(["  nlopt............................ : $enablenlopt"])
  AS_ECHO(["  parmetis......................... : $enableparmetis"])
  AS_ECHO(["  petsc............................ : $enablepetsc"])
  AS_IF([test "x$enablepetsc" = "xyes"],
        [echo '     'version....................... : $petscversion])
  AS_ECHO(["  qhull............................ : $enableqhull"])
  AS_ECHO(["  sfcurves......................... : $enablesfc"])
  AS_ECHO(["  slepc............................ : $enableslepc"])
  AS_IF([test "x$enableslepc" = "xyes"],
        [echo '     'version....................... : $slepcversion])
  AS_ECHO(["  thread model..................... : $found_thread_model"])
  AS_ECHO(["  c++ rtti ........................ : $ac_cv_cxx_rtti"])
  AS_ECHO(["  tecio............................ : $enabletecio"])
  AS_ECHO(["  tecplot...\(vendor binaries\)...... : $enabletecplot"])
  AS_ECHO(["  tetgen........................... : $enabletetgen"])
  AS_ECHO(["  triangle......................... : $enabletriangle"])
  AS_ECHO(["  trilinos......................... : $enabletrilinos"])
  AS_IF([test "x$enabletrilinos" = "xyes"],
        [
          AS_ECHO(["     AztecOO....................... : $enableaztecoo"])
          AS_ECHO(["     NOX........................... : $enablenox"])
          AS_ECHO(["     ML............................ : $enableml"])
          AS_ECHO(["     Tpetra........................ : $enabletpetra"])
          AS_ECHO(["     DTK........................... : $enabledtk"])
          AS_ECHO(["     Ifpack........................ : $enableifpack"])
          AS_ECHO(["     Epetra........................ : $enableepetra"])
          AS_ECHO(["     EpetraExt..................... : $enableepetraext"])
        ])
  AS_ECHO(["  vtk.............................. : $enablevtk"])
  AS_IF([test "x$enablevtk" = "xyes"],
        [echo '     'version....................... : $vtkversion])
  AS_ECHO([])
  AS_IF([test "x$libmesh_optional_INCLUDES" != "x"],
        [
          echo '  'libmesh_optional_INCLUDES........ : $libmesh_optional_INCLUDES
          AS_ECHO([])
        ])
  AS_IF([test "x$libmesh_optional_LIBS" != "x"],
        [
          echo '  'libmesh_optional_LIBS............ : $libmesh_optional_LIBS
          AS_ECHO([])
        ])
fi
echo '-------------------------------------------------------------------------------'

echo Configure complete, now type \'make\' and then \'make install\'.
AS_ECHO([])

])
