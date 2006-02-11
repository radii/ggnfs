<?php

$QPARAMS = "block_params";
$QFLAG = ".qflag";
$q0=7000000;
$qintsize=10000;
$qcount=1;
$addname="other";
$LOGFILENAME="alloc.log";

if( !array_key_exists( "name", $_POST) || !array_key_exists( "blocksize", $_POST ) || !array_key_exists( "blockcount", $_POST) )
{
 print("Incorrect parameters<br>\n");
 exit(0);
}

if( strlen($_POST["name"]) > 0 )
	$addname = $_POST["name"];

$addname = strtolower( $addname );
$addname = preg_replace("/\s+/", "", $addname);
$addname = preg_replace("/[^a-z0-9]+/", "", $addname);

$ip = $_SERVER["REMOTE_ADDR"];
$time = date("d.m.Y G:i:s");

$bsize = $_POST["blocksize"];
if ($bsize === "small") $qintsize=10000;
else if ($bsize === "medium") $qintsize=50000;
else if ($bsize === "large") $qintsize=250000;
else if ($bsize === "verylarge") $qintsize=500000;
else if ($bsize === "huge") $qintsize=1000000;

$qcount = intval($_POST["blockcount"]);

while( file_exists( $QFLAG) ) { usleep( 50000 ); };

touch($QFLAG);

$inf = fopen( $QPARAMS, "r");
chmod( $QPARAMS, 0777 );

$logf = fopen( $LOGFILENAME, "a");
chmod( $LOGFILENAME, 0777 );

$q0 = intval( fgets( $inf ) );
fclose($inf);

if( !file_exists( $addname ) )
{
 mkdir( $addname );
 chmod( $addname, 0777 );
}

for ($i=0; $i<$qcount; ++$i)
{
 print("<h2>Block # ".($i+1)."</h2><br>\n");
 $cq0 = $q0 + $i * $qintsize;
 fwrite( $logf, "$time,$ip,$addname,$cq0,$qintsize\n");
 $cdir="block-$cq0";
 mkdir("$addname/$cdir");
 chmod( "$addname/$cdir", 0777 );
 $outf = fopen( "$addname/$cdir"."/$cdir.job", "w");
 chmod( "$addname/$cdir"."/$cdir.job", 0777 );

fwrite( $outf, "name: distributed\n" );
fwrite( $outf, "type: gnfs\n" );
fwrite( $outf, "n: 68028769341463232121841916970039291388323866035762490102898304057696030390081374130780566873673063180687912398472130182736556136787752007\n" );
fwrite( $outf, "m: 19579601246936383255177337290943723213162756829393622310753821291194479478247833902142966211201826360563584005971061787712042388455520826\n" );
fwrite( $outf, "deg: 5\n");
fwrite( $outf, "skew: 126039.92\n" );
fwrite( $outf, "# norm 1.63e+19\n" );
fwrite( $outf, "c5: 5324400\n" );
fwrite( $outf, "c4: 2091687129165\n" );
fwrite( $outf, "c3: -155086943617402180\n" );
fwrite( $outf, "c2: -4778484126500310802610\n" );
fwrite( $outf, "c1: 485306265854886436022405128\n" );
fwrite( $outf, "c0: -38785275955005958598756653605031\n" );
fwrite( $outf, "# alpha -6.79\n" );
fwrite( $outf, "Y1: 1188447494830519\n" );
fwrite( $outf, "Y0: -105022898916743963791871982\n" );
fwrite( $outf, "# Murphy_E 3.01e-11\n" );
fwrite( $outf, "rlim: 11000000\n" );
fwrite( $outf, "alim: 11000000\n" );
fwrite( $outf, "lpbr: 27\n" );
fwrite( $outf, "lpba: 27\n" );
fwrite( $outf, "mfbr: 52\n" );
fwrite( $outf, "mfba: 52\n" );
fwrite( $outf, "rlambda: 2.5\n" );
fwrite( $outf, "alambda: 2.5\n" );
fwrite( $outf, "q0: $cq0\n");
fwrite( $outf, "qintsize: $qintsize\n");
fclose ($outf);

 print("<a href=\"http://factoring.org/distributed/$addname/$cdir/$cdir.job\" >$cdir.job</a><br>\n");

 $outf = fopen ( "$addname/$cdir"."/runme.cmd", "w" );
 chmod( "$addname/$cdir"."/runme.cmd", 0777 );
 fwrite( $outf, "start \"\" /LOW instr.set.cmd");
 fclose( $outf );

 print("<a href=\"http://factoring.org/distributed/$addname/$cdir/runme.cmd\" >runme.cmd</a><br>\n");

 $outf = fopen ( "$addname/$cdir"."/instr.set.cmd", "w" );
 chmod( "$addname/$cdir"."/instr.set.cmd", 0777 );
 fwrite( $outf, "call \"gnfs-lasieve4I14e.exe\" -k -o spairs.$cq0-".($cq0 + $qintsize).".out -v -a block-$cq0.job -n1234\n");
 fwrite( $outf, "call tar -cf result-$cq0.tar spairs.$cq0-".($cq0 + $qintsize).".out ggnfs.log\n");
 fwrite( $outf, "call bzip2 -9 result-$cq0.tar\n");
 fclose ($outf);

 print("<a href=\"http://factoring.org/distributed/$addname/$cdir/instr.set.cmd\" >instr.set.cmd</a><br>\n");

}

$inf = fopen ( $QPARAMS, "w+" );

$q0 = $q0 + $qintsize * $qcount;

fwrite( $inf, "$q0");
fclose ($inf);
fclose ($logf);

unlink($QFLAG);

?>