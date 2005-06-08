l1_bits		.ASSIGNA	15
longer_ri	.ASSIGNA	0
lasched_updates	.ASSIGNA	1
n_i		.ASSIGNA	1
.AREPEAT	\&i_bits
n_i		.ASSIGNA	2 * \&n_i
.AENDR
n_i_mask	.ASSIGNA	\&n_i - 1
l1_size		.ASSIGNA	1
.AREPEAT	\&l1_bits
l1_size		.ASSIGNA	2 * \&l1_size
.AENDR
l1_mask		.ASSIGNA	\&l1_size - 1
j_per_strip	.ASSIGNA	\&l1_size / \&n_i
