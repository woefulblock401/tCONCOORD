## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.


AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\" 


LDADD = libcnc@LIBSUFFIX@.la ../mdlib/libmd@LIBSUFFIX@.la \
../gmxlib/libgmx@LIBSUFFIX@.la 


lib_LTLIBRARIES = libcnc@LIBSUFFIX@.la


libcnc@LIBSUFFIX@_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@


libcnc@LIBSUFFIX@_la_SOURCES = \
	tconcoord.h  disco_util.c  atomlist.c  resl.c histogram.c  \
	cnctop.c    array.c     contab.c     bounds.c  boundtrack.c lib_util.c \
	  dist_util.c  planar.c      error.c      chiral.c    dihed.c \
  structs.h vdwcomb.c   hbonds.c      nbsearch.c   nonbonded.c  \
  distance.c   bonded.c    dihedral.c  structure.c  vdw.c      \
  groups.c     namegroups.c  geometry.c   bondlib.c   \
  angle.c     string_util.c  random.c        types.c  \
	input.c       force.c      bondlist.c  exfo.c sidechain.c \
rosetta.c pdbqt.c gridmap.c packing.c buried.c rotation.c sasa.c rama.c

bin_PROGRAMS = \
	tdist tdisco tsearch tinteract  packscore doplanar hb_anal fuse_xtc bsfind tprog solv_energy score_hbonds



tdisco_SOURCES = tdisco.c


tdist_SOURCES = tdist.c

tsearch_SOURCES = tsearch.c

bsfind_SOURCES = bsfind.c

tinteract_SOURCES  = tinteract.c

packscore_SOURCES = packscore.c

doplanar_SOURCES = doplanar.c

hb_anal_SOURCES = hb_anal.c

fuse_xtc_SOURCES = fuse_xtc.c

tprog_SOURCES = tprog.c

score_hbonds_SOURCES = score_hbonds.c

CLEANFILES = *~ \\\#*


