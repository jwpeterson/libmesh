BUILT_SOURCES = exact_solution.C
DISTCLEANFILES = $(BUILT_SOURCES)
exact_solution.C: $(top_srcdir)/examples/introduction/introduction_ex3/exact_solution.C
	$(AM_V_GEN)$(LN_S) -f $< $@

example_name  = subdomains_ex2
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/subdomains/ex2
data          = subdomains_ex2.C exact_solution.C run.sh
sources       = $(data)

CLEANFILES = out_*.e gnuplot_*

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
