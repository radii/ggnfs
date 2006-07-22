@echo off
rem **************************************************************************
rem * relsort.cpp
rem * Copyleft 2006 by Sten
rem *
rem * This is batch file to automate relsort.exe build process
rem * on Win32 platform under MinGW. You should have MinGW and MSYS
rem * installed.
rem * Typical build sequence looks like:
rem *
rem *  make config
rem *  make all
rem *  make clean
rem *
rem * It is free software; you can redistribute it and/or modify
rem * it under the terms of the GNU General Public License as published by
rem * the Free Software Foundation; either version 2 of the License, or
rem * (at your option) any later version.
rem *
rem *************************************************************************

if "%1"=="" goto usage
if "%1"=="/?" goto usage
if "%1"=="-?" goto usage
if "%1"=="\?" goto usage
if "%1"=="-help" goto usage
if "%1"=="/help" goto usage
if "%1"=="all" goto all
if "%1"=="config" goto config
if "%1"=="clean" goto clean

:usage
echo.
echo   Usage: "make all|config|clean]"
echo.
exit

:config
set SAVEDDIR=%CD%
cd ../glob
sh configure 
cd %SAVEDDIR%
exit;

:all
set SAVEDDIR=%CD%
cd ../glob
make 
gcc -c -O3 -DNDEBUG fnmatch.c 
cd %SAVEDDIR%
g++ -s -O3 -DNDEBUG -I../glob relsort.cpp ../glob/fnmatch.o ../glob/glob.o -o relsort.exe
exit;

:clean
del *.o
del config.log
del relsort.exe
set SAVEDDIR=%CD%
cd ../glob
make clean
cd %SAVEDDIR%
exit;
