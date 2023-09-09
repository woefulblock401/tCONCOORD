#!/bin/csh

setenv root `pwd`
setenv srcdir $root/src/tconcoord
mkdir tcnc_bin_dist
cd tcnc_bin_dist
mkdir bin

foreach S(`cat $root/tcnc.binarylist`)
echo Copying $srcdir/$S
cp $srcdir/$S bin/
end

cp -r $root/cnclib  lib
cp -r $root/tcnc_examples examples

foreach S(`find . -name ".svn"`)
rm -rf $S
end
cd ..

