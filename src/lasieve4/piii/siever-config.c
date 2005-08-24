/* siever-config.c
  6/13/04: Hacked up for inclusion in GGNFS by Chris Monico.
                                                                                                       
  Copyright (C) 2001 Jens Franke, T. Kleinjung.
  This file is part of gnfs4linux, distributed under the terms of the
  GNU General Public Licence and WITHOUT ANY WARRANTY.
                                                                                                       
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
  */
#include "siever-config.h"
                                                                                                       
void siever_init(void)
{
  init_montgomery_multiplication();
}
const uint32_t schedule_primebounds[N_PRIMEBOUNDS] =
  { 0x100000, 0x200000, 0x400000, 0x800000, 0x1000000, 0x2000000, ULONG_MAX };

const uint32_t schedule_sizebits[N_PRIMEBOUNDS] = { 20, 21, 22, 23, 24, 25, 32 };
