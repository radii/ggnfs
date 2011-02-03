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
dnl
define(aux1,r14)dnl
define(aux2,r15)dnl
define(aux3,r16)dnl
define(xc,r8)dnl
define(yc,r9)dnl
define(x,r32)dnl
define(y,r33)dnl
define(A,r34)dnl
define(auxf1,f6)dnl
define(auxf2,f7)dnl
define(auxf3,f8)dnl
define(auxf4,f9)dnl
define(a1,f10)dnl
define(a2,f11)dnl
define(quot,f12)dnl
define(q1,f13)dnl
define(corr,f14)dnl
define(eps,f15)dnl
define(rec,f127)dnl
define(nstep,5)dnl
dnl
define(euclid_step,
	add $2=$2`,'$4
	sub $1=$1`,'$3;;
	cmp.ltu p6`,'p7=$1`,'$3;;
	(p7)add $2=$2`,'$4
	(p7)sub $1=$1`,'$3;;
	cmp.ltu p6`,'p7=$1`,'$3;;
	(p7)add $2=$2`,'$4
	(p7)sub $1=$1`,'$3;;
	cmp.ltu p6`,'p7=$1`,'$3;;
	(p7)add $2=$2`,'$4
	(p7)sub $1=$1`,'$3;;
	cmp.ltu p6`,'p7=$1`,'$3;;
	(p6)br.sptk ``es_ende''es_count
	add $2=$2`,'$4
	sub $1=$1`,'$3
	shl aux1=$3`,'nstep
	shl aux2=$4`,'nstep
	shl aux3=$3`,'eval(nstep+1);;
	cmp.ltu p6`,'p7=$1`,'aux3
	cmp.leu p8`,'p9=aux1`,'$1;;
	(p7)br.spnt ``float_euclid''es_count
	forloop(`i',0,nstep,
		(p8)sub $1=$1``,''aux1
		shr.u aux1=aux1``,''1
		(p8)add $2=$2``,''aux2 ;;
		cmp.leu p8``,''p9=aux1``,''$1
		shr.u aux2=aux2``,''1;;
	)
	br.sptk.many ``es_ende''es_count
``float_euclid''es_count:
	setf.sig auxf3=$3
	setf.sig auxf1=$1;;
	setf.sig auxf2=$2
	fcvt.xf auxf3=auxf3
	fcvt.xf auxf1=auxf1;;
	mov a1=auxf1
	setf.sig auxf4=$4
	fcvt.xf auxf2=auxf2;;
	fcvt.xf auxf4=auxf4
	frcpa.s1 rec`,'p6=auxf1`,'auxf3;;
	fma.s1	quot=auxf1`,'rec`,'f0
	fnma.s1 corr=rec`,'auxf3`,'f1;;
	fcvt.fxu.trunc.s1 q1=quot
	fma.s1 quot=corr`,'quot`,'quot;;
	fma.s1 corr=corr`,'corr`,'eps
	fcvt.xf q1=q1;;
	fnma.s1 auxf1=q1`,'auxf3`,'auxf1
	mov a2=auxf2;;
	fma.s1 auxf2=q1`,'auxf4`,'auxf2
	fcmp.le p6`,'p7=f0`,'auxf1;;
	(p6)fcmp.lt.unc p6`,'p7=auxf1`,'auxf3
	fma.s1 quot=corr`,'quot`,'quot;;
	(p6)br.sptk.many ``fes_ende''es_count
	fcvt.fxu.trunc.s1 q1=quot;;
	fcvt.xf q1=q1;;
	fnma.s1 auxf1=q1`,'auxf3`,'a1
	fma.s1 auxf2=q1`,'auxf4`,'a2;;
``fes_ende''es_count:
	fcvt.fxu.s1 auxf1=auxf1
	fcvt.fxu.s1 auxf2=auxf2;;
	getf.sig $1=auxf1
	getf.sig $2=auxf2;;
``es_ende''es_count:`define(`es_count',eval(es_count+1))')dnl
define(es_count,2)dnl esteps 0 and 1 in modinv, 2-5 in ri-aux
function_head(asm_getbc)
	alloc r31=ar.pfs,7,0,0,0
	zxt4 x=x
	mov r2=0xffdd;;
	zxt4 y=y
	cmp.ltu p6,p7=x,A;;
 	add xc=1,r0
	mov yc=r0
	(p6)br.dpnt have_bs
	setf.exp eps=r2
divide:
	euclid_step(y,yc,x,xc)
	cmp.ltu p6,p7=y,A;;
	(p6)br.dpnt have_ct
	euclid_step(x,xc,y,yc)
	cmp.ltu p6,p7=x,A;;
	(p7)br.dptk divide
have_bs:
	# Have (x,xc)=(b,s). One final Euclid step to get (c,t)
	st4 [r35]=x
	sub A=A,x;;
	st4 [r36]=xc
	sub y=y,A;;
	euclid_step(y,yc,x,xc)
	st4 [r38]=yc
	add y=y,A;;
	st4 [r37]=y
	br.ret.sptk b0
have_ct:
	# Have (y,yc)=(c,t). One final Euclid step to get (b,s)
	st4 [r37]=y
	sub A=A,y;;
	st4 [r38]=yc
	sub x=x,A;;
	euclid_step(x,xc,y,yc)
	st4 [r36]=xc
	add x=x,A;;
	st4 [r35]=x
	br.ret.sptk b0
