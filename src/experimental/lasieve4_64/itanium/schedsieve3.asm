define(function_head,
.text
	.align 16
	.globl $1#
	.proc $1#
$1:
)dnl
dnl This macro is taken from the GNU info documentation of m4.
define(`forloop',
       `pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')dnl
define(`_forloop',
       `$4`'ifelse($1, `$3', ,
     		 `define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')dnl
define(slog,r32)dnl
define(sieve,r33)dnl
define(sched,r34)dnl
define(sched_ub,r35)dnl
define(advanced_sched,r8)dnl
define(sv0,r16)dnl
define(sv1,r17)dnl
define(sv2,r18)dnl
define(sv3,r19)dnl
define(sv4,r20)dnl
define(sv5,r21)dnl
define(sv6,r22)dnl
define(sv7,r23)dnl
define(sloc0,r24)dnl
define(sloc1,r25)dnl
define(sloc2,r26)dnl
define(sloc3,r27)dnl
define(sloc4,r28)dnl
define(sloc5,r29)dnl
define(sloc6,r30)dnl
define(sloc7,r31)dnl
function_head(schedsieve)
	alloc r31=ar.pfs,4,0,0,0
	add sched_ub=-28,sched_ub;;
	cmp.ltu p6,p7=sched,sched_ub
	add advanced_sched=288,sched
	ld2 sloc0=[sched],4;;
	ld2 sloc1=[sched],4;;
	ld2 sloc2=[sched],4;;
	ld2 sloc3=[sched],4;;
	ld2 sloc4=[sched],4;;
	add sloc0=sloc0,sieve
	ld2 sloc5=[sched],4;;
	ld1.a sv0=[sloc0]
	add sloc1=sloc1,sieve
	ld2 sloc6=[sched],4;;
	ld2 sloc7=[sched],4
	add sloc2=sloc2,sieve
	ld1.a sv1=[sloc1]
	add sv0=sv0,slog
	(p7)br.spnt schedsieve_loopende;;
schedsieve_loop:
	cmp.ltu p6,p7=sched,sched_ub
	lfetch [advanced_sched],32

	ld1.a sv2=[sloc2]
	add sloc3=sloc3,sieve
	chk.a.clr sv0,reread_sloc0
have_sv0:
	st1 [sloc0]=sv0
	add sv1=sv1,slog
	ld2 sloc0=[sched],4;;

	ld1.a sv3=[sloc3]
	add sloc4=sloc4,sieve
	chk.a.clr sv1,reread_sloc1
have_sv1:
	st1 [sloc1]=sv1
	add sv2=sv2,slog
	ld2 sloc1=[sched],4;;

	ld1.a sv4=[sloc4]
	add sloc5=sloc5,sieve
	chk.a.clr sv2,reread_sloc2
have_sv2:
	st1 [sloc2]=sv2
	add sv3=sv3,slog
	ld2 sloc2=[sched],4;;

	ld1.a sv5=[sloc5]
	add sloc6=sloc6,sieve
	chk.a.clr sv3,reread_sloc3
have_sv3:
	st1 [sloc3]=sv3
	add sv4=sv4,slog
	ld2 sloc3=[sched],4;;

	ld1.a sv6=[sloc6]
	add sloc7=sloc7,sieve
	chk.a.clr sv4,reread_sloc4
have_sv4:
	st1 [sloc4]=sv4
	add sv5=sv5,slog
	ld2 sloc4=[sched],4;;

	ld1.a sv7=[sloc7]
	add sloc0=sloc0,sieve
	chk.a.clr sv5,reread_sloc5
have_sv5:
	st1 [sloc5]=sv5
	add sv6=sv6,slog
	ld2 sloc5=[sched],4;;

	ld1.a sv0=[sloc0]
	add sloc1=sloc1,sieve
	chk.a.clr sv6,reread_sloc6
have_sv6:
	st1 [sloc6]=sv6
	add sv7=sv7,slog
	ld2 sloc6=[sched],4;;

	ld1.a sv1=[sloc1]
	add sloc2=sloc2,sieve
	chk.a.clr sv7,reread_sloc7
have_sv7:
	st1 [sloc7]=sv7
	add sv0=sv0,slog
	ld2 sloc7=[sched],4;;

	(p6)br.sptk schedsieve_loop
schedsieve_loopende:
	add sloc3=sloc3,sieve
	add sloc4=sloc4,sieve
	add sloc5=sloc5,sieve
	add sloc6=sloc6,sieve
	add sloc7=sloc7,sieve
	add sched_ub=28,sched_ub
	add sched=-32,sched;;
forloop(`i',0,6,`define(`xyz',
	cmp.ltu p6`,'p7=sched`,'sched_ub;;
	(p6)ld1 sv0=[sloc`'i]
	add sched=4`,'sched;;
	(p6)add sv0=sv0`,'slog;;
	(p6)st1 [sloc`'i]=sv0;;
)xyz')
	br.ret.sptk b0
forloop(`i',0,7,`define(`xyz',
	.align 16
reread_sloc`'i:
	ld1 sv`'i=[sloc`'i];;
	add sv`'i=sv`'i`,'slog
	br.sptk have_sv`'i;;
)xyz')
