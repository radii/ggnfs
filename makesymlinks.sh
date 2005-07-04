#!/bin/sh

ln -s src/Makefile.athlon src/Makefile
ln -s srs/lasieve4/piii src/lasieve4/asm
chmod a+x src/autogplot.sh
chmod a+x tests/factLat.pl
echo "Symbolic links have been created and execute permissions have been fixed."
