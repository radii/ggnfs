	.globl schedsieve
/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/

schedsieve:
	pushl %esi
	pushl %edi
	pushl %ebx
	movl 28(%esp),%edx
	movb 16(%esp),%al
	movl 20(%esp),%edi
	subl $12,%edx
	movl 24(%esp),%esi
	cmpl %esi,%edx
        movzwl (%esi),%ebx
        jbe fat_loop_end
fat_loop:
#if 1
        prefetcht1 128(%esi)
#endif
        movzwl 4(%esi),%ecx
        addb %al,(%edi,%ebx)
        movzwl 8(%esi),%ebx
        addb %al,(%edi,%ecx)
        movzwl 12(%esi),%ecx
        leal 16(%esi),%esi
        addb %al,(%edi,%ebx)
        movzwl (%esi),%ebx
        addb %al,(%edi,%ecx)
        cmpl %esi,%edx
        ja fat_loop
fat_loop_end:
	addl $12,%edx
	cmpl %esi,%edx
        movzwl 4(%esi),%ecx
	leal 4(%esi),%esi
	jbe schedsieve_end
        addb %al,(%edi,%ebx)
	cmpl %esi,%edx
        movzwl 4(%esi),%ebx
	leal 4(%esi),%esi
	jbe schedsieve_end
        addb %al,(%edi,%ecx)
	cmpl %esi,%edx
	jbe schedsieve_end
        addb %al,(%edi,%ebx)
schedsieve_end:
	popl %ebx
	popl %edi
	popl %esi
	ret
