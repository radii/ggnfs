define(l1_bits,16)dnl
define(n_i,eval(2**n_i_bits))dnl
define(n_i_mask,eval(n_i-1))dnl
define(l1_size,eval(2**l1_bits))dnl
define(l1_mask,eval(l1_size-1))dnl
define(j_per_strip,eval(2**(l1_bits-n_i_bits)))dnl
