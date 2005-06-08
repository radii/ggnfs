#!/bin/bash
USAGE="USAGE: $0 [jpg output name] [ <data1 name> <data1 file> ] ... "
TMPFILE="tmpfile.gp"
GP=/usr/bin/gnuplot

if [ ! -x $GP ] ; then
  echo "$0 : Cannot find gnuplot. Terminating..."
  exit 0
fi

typeset -i NAMELOC=2
typeset -i DATALOC=3
typeset -i TMPINT=0

if test -z $3 ; then
  echo $USAGE
  exit
fi
if test -e $TMPFILE; then
  echo "Tempfile $TMPFILE already exists! Cannot continue!"
  exit
fi
#echo "set terminal jpeg transparent interlace small size 1024,800" >> $TMPFILE
echo "set terminal png small " >> $TMPFILE
echo "set out '$1'" >> $TMPFILE
if test -n '$XAXIS'; then
  echo "set xlabel '$XAXIS'" >> $TMPFILE
else  
  echo "set xlabel 'Total relations'" >> $TMPFILE
fi;
if test -n '$YAXIS'; then
  echo "set ylabel '$YAXIS'" >> $TMPFILE
else
  echo "set ylabel 'Full relations'" >> $TMPFILE
fi

#echo  "set format y \"10^{%L}\"" >> $TMPFILE
#echo  "set logscale y" >> $TMPFILE

PRE="plot "
COMMAND=""
while test $DATALOC -le $# ; do
  TRAILER=""
  TMPINT=$DATALOC+2
  if test $TMPINT -le $# ; then
    TRAILER=", " ;
  fi ;
  COMMAND="$COMMAND $PRE '${!DATALOC}' using 1:2 title '${!NAMELOC}' with lines $TRAILER "
  PRE=" "
  DATALOC=($DATALOC+2)
  NAMELOC=($NAMELOC+2)
done
echo $COMMAND >> $TMPFILE
$GP $TMPFILE
rm -f $TMPFILE

