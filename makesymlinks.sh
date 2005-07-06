#!/bin/sh

cd src; ln -s Makefile.athlon Makefile; cd ..
cd src/lasieve4; ln -s piii asm; cd ../..
chmod a+x src/autogplot.sh
chmod a+x tests/factLat.pl
echo "Symbolic links have been created and execute permissions have been fixed."
