# csh/tcsh configuration file for Gromacs.
# First we remove previous Gromacs stuff from paths 
# by selecting everything else. 
# Removal is not 100% necessary, but very useful when we
# repeatedly switch between gmx versions in a shell.

# zero possibly unset vars to avoid warnings
if (! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH ""
if (! $?PATH) setenv PATH ""
if (! $?MANPATH) setenv MANPATH ""
if (! $?GMXLDLIB) setenv GMXLDLIB ""
if (! $?GMXBIN) setenv GMXBIN ""
if (! $?GMXMAN) setenv GMXMAN ""

# remove previous gromacs part from ld_library_path
set tmppath = ""
foreach i ( `echo $LD_LIBRARY_PATH | sed "s/:/ /g"` )
  if ( "$i" != "$GMXLDLIB" ) set tmppath = "${tmppath}:$i"
end
setenv LD_LIBRARY_PATH $tmppath

# remove gromacs stuff from binary path
set tmppath = ""
foreach i ( `echo $PATH | sed "s/:/ /g"` )
  if ( "$i" != "$GMXBIN" ) set tmppath = "${tmppath}:$i"
end
setenv PATH $tmppath

# and remove stuff from manual path
set tmppath = ""
foreach i ( `echo $MANPATH | sed "s/:/ /g"` )
  if ( "$i" != "$GMXMAN" ) set tmppath = "${tmppath}:$i"
end
setenv MANPATH $tmppath

##########################################################
# This is the real configuration part. We save the Gromacs
# things in separate vars, so we can remove them later.
# If you move gromacs, change the next four line.
##########################################################
setenv GMXBIN /usr/local/gromacs/bin
setenv GMXLDLIB /usr/local/gromacs/lib
setenv GMXMAN /usr/local/gromacs/share/man
setenv GMXDATA /usr/local/gromacs/share

# old variables begin with ':' now, or are empty.
setenv PATH ${GMXBIN}${PATH}
setenv LD_LIBRARY_PATH ${GMXLDLIB}${LD_LIBRARY_PATH}
setenv MANPATH ${GMXMAN}${MANPATH}

setenv GMXFONT	10x20

# Read completions if we understand it (i.e. have tcsh)
if { complete >& /dev/null } then
  if ( -f $GMXBIN/completion.csh ) source $GMXBIN/completion.csh
endif








