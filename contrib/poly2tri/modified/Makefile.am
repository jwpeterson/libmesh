# [JWP] pkg_sources and pkg_libs are not automake-specific variables as far as I can tell...
# These variables will thus go through unmodified to the Makefile generated from this Makefile.am file.
pkg_sources = $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/shapes.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/advancing_front.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/cdt.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep_context.cc \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/poly2tri.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/shapes.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/utils.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/advancing_front.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/cdt.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep.h \
              $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep_context.h
pkg_libs    =

# In Automake, EXTRA_DIST is a variable used to specify files or
# directories that should be included in the source distribution
# created by make dist, even if they are not explicitly referenced by
# other Automake rules or automatically detected.
EXTRA_DIST = $(top_srcdir)/contrib/poly2tri/poly2tri/AUTHORS \
             $(top_srcdir)/contrib/poly2tri/poly2tri/LICENSE \
             $(top_srcdir)/contrib/poly2tri/poly2tri/README.md

# import as much commonality as possible
include $(top_srcdir)/contrib/Make.common

# These are the headers that we "modify" from the original source code
# by prefixing the #include statments with the "poly2tri" path and
# copying them into the "contrib/poly2tri/modified" sub-directory of
# the build directory. This is not an variable recognized by automake,
# see BUILT_SOURCES below for that.
poly2tri_hfiles = poly2tri/poly2tri.h \
                  poly2tri/common/shapes.h \
                  poly2tri/common/utils.h \
                  poly2tri/sweep/advancing_front.h \
                  poly2tri/sweep/cdt.h \
                  poly2tri/sweep/sweep_context.h \
                  poly2tri/sweep/sweep.h

# This is an automake variable which defines a list of header files
# that are built and must not be distributed. [JWP] I think this refers
# to files that should not be included in tar archives created by
# "make dist" and is not related to installation.
nodist_include_HEADERS = $(poly2tri_hfiles)

# This is an automake variable that tells what headers are installed.
# The "nobase" prefix maintains the relative directory structure of
# the headers when they are installed.
# https://www.gnu.org/software//automake/manual/html_node/Headers.html
nobase_include_HEADERS = $(poly2tri_hfiles)

# This is an automake variable that lists headers used by programs or
# convenience libraries that should not be installed.
noinst_HEADERS = $(poly2tri_hfiles)

# BUILD_SOURCES is an automake variable that is used to specify source
# files that need to be generated early in the build process but are
# not compiled themselves into the final executable or library
BUILT_SOURCES = $(poly2tri_hfiles)

clean-local:
	rm -f $(BUILT_SOURCES)

poly2tri/poly2tri.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/poly2tri.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/
	$(AM_V_GEN) $(SED) -e 's#include "#include "poly2tri/#' < $< > $@

poly2tri/common/shapes.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/shapes.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/common/
	$(AM_V_GEN) rm -f $@ && cp -f $< $@

poly2tri/common/utils.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/common/utils.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/common/
	$(AM_V_GEN) $(SED) -e 's#include "#include "poly2tri/#' < $< > $@

poly2tri/sweep/advancing_front.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/advancing_front.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/sweep/
	$(AM_V_GEN) $(SED) -e 's#include "../#include "poly2tri/#' < $< > $@

poly2tri/sweep/cdt.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/cdt.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/sweep/
	$(AM_V_GEN) $(SED) -e 's#include "#include "poly2tri/sweep/#' < $< > $@

poly2tri/sweep/sweep_context.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep_context.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/sweep/
	$(AM_V_GEN) rm -f $@ && cp -f $< $@

poly2tri/sweep/sweep.h : $(top_srcdir)/contrib/poly2tri/poly2tri/poly2tri/sweep/sweep.h
	$(AM_V_GEN) $(MKDIR_P) poly2tri/sweep/
	$(AM_V_GEN) rm -f $@ && cp -f $< $@
