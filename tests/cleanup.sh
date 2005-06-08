#!/bin/sh
# This is just for cleaning up subdirectores prior to tarballing.
for file in * ; do
  if test -d $file ; then
    echo "Cleaning directory $file..."
    make -C $file squeaky > /dev/null 2>&1
  fi
done

