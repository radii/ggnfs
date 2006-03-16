if not exist %1 goto exit
if not exist %2 goto copy
fc %1 %2 > nul
if not %errorlevel 1 goto exit
:copy
echo copying %1 to %2 
copy %1 %2
:exit
