#!/bin/sh

../pol51m0b -b rsa100 -p 4 -a 1 -A 50 -n 5e14 -v
dos2unix rsa100.51.m
diff rsa100.51.m.old rsa100.51.m
rm rsa100.51.m
