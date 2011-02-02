

define(forloop,  
pushdef($1, $2)_forloop($1, $2, $3, $4)popdef($1))  
define(_forloop,  
$4@ifelse($1, $3, ,  
define($1, incr($1))_forloop($1, $2, $3, $4)))  
define(function_head,.text  
.p2align 4,,15  
.globl $1  
.type	$1, @function  
$1:)dnl  