#!/bin/sh

../pol51m0b -b rsa100 -p 4 -a 1 -A 50 -n 5e14 -v
dos2unix rsa100.51.m
diff rsa100.51.m.old rsa100.51.m
rm rsa100.51.m

../pol51m0b -b rsa110 -p 4 -a 1 -A 10 -n 2.52e16 -v
dos2unix rsa110.51.m
diff rsa110.51.m.old rsa110.51.m
rm rsa110.51.m

../pol51m0b -b rsa120 -p 4 -a 9 -A 10 -n 1.14e18 -v
dos2unix rsa120.51.m
diff rsa120.51.m.old rsa120.51.m
rm rsa120.51.m

../pol51m0b -b rsa130 -p 5 -a 27 -A 30 -n 5.16e19 -v
dos2unix rsa130.51.m
diff rsa130.51.m.old rsa130.51.m
rm rsa130.51.m

../pol51m0b -b rsa140 -p 6 -a 44 -A 50 -n 2.34e21 -v
dos2unix rsa140.51.m
diff rsa140.51.m.old rsa140.51.m
rm rsa140.51.m

