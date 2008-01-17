#!/bin/sh

../pol51opt -b rsa100_2 -n 2.86e13 -N 3.56e11 -e 2.8e-9 -v
dos2unix rsa100_2.cand
diff rsa100_2.cand.old rsa100_2.cand
rm rsa100_2.cand
