example_name  = vector_fe_ex3
check_SCRIPTS = run.sh # Convergence too poor with default solvers
install_dir   = $(examples_install_path)/vector_fe/ex3
data          = curl_curl_exact_solution.h vector_fe_ex3.in curl_curl_system.C solution_function.h curl_curl_system.h vector_fe_ex3.C run.sh
sources       = $(data)

CLEANFILES = out.e

# also need to link files for VPATH builds
if LIBMESH_VPATH_BUILD
  BUILT_SOURCES = .linkstamp
.linkstamp:
	-rm -f vector_fe_ex3.in && $(LN_S) -f $(srcdir)/vector_fe_ex3.in .
	$(AM_V_GEN)touch .linkstamp

  CLEANFILES += vector_fe_ex3.in .linkstamp
endif

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
