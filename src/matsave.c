/**************************************************************
 * matsave.c
 * Save matsolve state at regular intervals and when interrupted.
 * Copyright 2006 Geoffrey Reynolds
 **************************************************************/

/*  This file is part of GGNFS.
 *
 *   GGNFS is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   GGNFS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with GGNFS; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* TODO:
    Make stronger checks that the save file matches the loaded matrix.
    Save the Lanczos seed in the save file, and change matsolve.c so
     that the choice of seed is not reported until after resuming.
*/

#include <signal.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include "ggnfs.h"
#include "if.h"

/* /etc/magic entry:

   0       ulong   0x4040513B      GGNFS matsolve save file
   >4      ulong   <3              version %lu,
   >>8     long    x               %ld columns,
   >>12    ulong   x               iteration %lu.
   >4      ulong   >2              version %lu.
*/
#define MATSAVE_MAGIC 0x4040513B
#define MATSAVE_VERSION 2
#define MATSAVE_FILE_NAME "matsave"
#define MATSAVE_BACKUP_NAME "matsave.bak"

/* Exported variables */
volatile int matsave_interval = 300; /* By default save each 300 seconds */

/* Local variables */
static u32 matsave_magic = MATSAVE_MAGIC;
static u32 matsave_version = MATSAVE_VERSION;
static FILE *file = NULL;
static int using_backup = 0;

static void handle_signal(int signum)
{
  switch (signum) {
    case SIGTERM:
    case SIGINT:
      matsave_interval = -1;
      break;
  }
}

/* Set up signal handlers. */
static void init_matsave(void)
{
  signal(SIGINT, handle_signal);
  signal(SIGTERM, handle_signal);
}

static const char *file_name(void)
{
  return using_backup ? MATSAVE_BACKUP_NAME : MATSAVE_FILE_NAME;
}

static const char *file_desc(void)
{
  return using_backup ?
    "backup save file '" MATSAVE_BACKUP_NAME "'" :
    "save file '" MATSAVE_FILE_NAME "'";
}

static int errorclose0(const char *msg)
{
  fprintf(stderr, "%s %s.\n", msg, file_desc());
  fclose(file);
  return 0;
}

static int matread_u32(u32 *dst, size_t count)
{
  if (read_u32(file, dst, count) != count)
    return errorclose0("Error while reading from");
  return 1;
}

static int matread_s32(s32 *dst, size_t count)
{
  if (read_i32(file, dst, count) != count)
    return errorclose0("Error while reading from");
  return 1;
}

static int matread_u64(u64 *dst, size_t count)
{
  if (read_u64(file, dst, count) != count)
    return errorclose0("Error while reading from");
  return 1;
}

static int matwrite_u32(const u32 *src, size_t count)
{
  if (write_u32(file, src, count) != count)
    return errorclose0("Error while writing to");
  return 1;
}

static int matwrite_s32(const s32 *src, size_t count)
{
  if (write_i32(file, src, count) != count)
    return errorclose0("Error while writing to");
  return 1;
}

static int matwrite_u64(const u64 *src, size_t count)
{
  if (write_u64(file, src, count) != count)
    return errorclose0("Error while writing to");
  return 1;
}

/* Attempt to resume from the save file or, failing that, from the
   backup. Return the current iteration if successful, 0 if not. */
u32 matresume(s32 n, u64 *Wi, u64 *Wi_1, u64 *Wi_2, u64 *T_1,
              u64 *tmp, u64 *U_1, u64 *tmp2, int *Si, int *Si_1,
              u64 *X, u64 *Y, u64 *Vi, u64 *Vi_1, u64 *Vi_2)
{
  u32 magic, version, iterations;
  s32 columns;

  init_matsave();
  for (using_backup = 0; using_backup < 2; using_backup++) {
    if ((file = fopen(file_name(),"rb")) == NULL) {
      if (errno != ENOENT) /* missing save file is not an error. */
        fprintf(stderr, "Could not open %s for reading.\n", file_desc());
      continue;
    }
    if (!matread_u32(&magic,1)) continue;
    if (magic != matsave_magic) {
      errorclose0("Wrong magic number in");
      continue;
    }
    if (!matread_u32(&version,1)) continue;
    if (version != matsave_version) {
      fprintf (stderr, "Could not read version %" PRIu32 " %s.\n",
               version, file_desc());
      fclose(file);
      continue;
    }
    if (!matread_s32(&columns,1)) continue;
    if (n != columns) {
      errorclose0("Loaded matrix does not match that in");
      continue;
    }
    if (!matread_u32(&iterations,1)) continue;
    if (!matread_u64(Wi,64)) continue;
    if (!matread_u64(Wi_1,64)) continue;
    if (!matread_u64(Wi_2,64)) continue;
    if (!matread_u64(T_1,64)) continue;
    if (!matread_u64(tmp,64)) continue;
    if (!matread_u64(U_1,64)) continue;
    if (!matread_u64(tmp2,64)) continue;
    /* Si and Si_1 are arrays of int, but stored as arrays of s32. */
    if (sizeof(int) != sizeof(s32)) {
      s32 i, sbuf[64];
      if (!matread_s32(sbuf,64)) continue;
      for (i = 0; i < 64; i++)
        Si[i] = sbuf[i];
      if (!matread_s32(sbuf,64)) continue;
      for (i = 0; i < 64; i++)
        Si_1[i] = sbuf[i];
    }
    else {
      if (!matread_s32(Si,64)) continue;
      if (!matread_s32(Si_1,64)) continue;
    }
    if (!matread_u64(X,n)) continue;
    if (!matread_u64(Y,n)) continue;
    if (!matread_u64(Vi,n)) continue;
    if (!matread_u64(Vi_1,n)) continue;
    if (!matread_u64(Vi_2,n)) continue;
    fclose(file);
    fprintf(stdout, "Resuming from %s at iteration %" PRIu32 ".\n",
            file_desc(), iterations);
    return iterations;
  }
  return 0;
}

/* Attempt to create a save file. Return 1 if successful, 0 if not. */
int matsave(u32 iterations, s32 n, const u64 *Wi, const u64 *Wi_1,
            const u64 *Wi_2, const u64 *T_1, const u64 *tmp,
            const u64 *U_1, const u64 *tmp2, const int *Si,
            const int *Si_1, const u64 *X, const u64 *Y,
            const u64 *Vi, const u64 *Vi_1, const u64 *Vi_2)
{
  using_backup = 0;
  /* Sten: delete matsave.bak before renaming matsave to matsave.bak. */
  remove(MATSAVE_BACKUP_NAME); 

  if (rename(MATSAVE_FILE_NAME,MATSAVE_BACKUP_NAME) == -1 && errno != ENOENT)
    fprintf(stderr, "Could not backup %s.\n", file_desc());
  if ((file = fopen(MATSAVE_FILE_NAME,"wb")) == NULL) {
    fprintf(stderr, "Could not open %s for writing.\n", file_desc());
    return 0;
  }

  printTmp("Creating %s at iteration %" PRIu32 "...",file_desc(),iterations);
  if (!matwrite_u32(&matsave_magic,1)) return 0;
  if (!matwrite_u32(&matsave_version,1)) return 0;
  if (!matwrite_s32(&n,1)) return 0;
  if (!matwrite_u32(&iterations,1)) return 0;
  if (!matwrite_u64(Wi,64)) return 0;
  if (!matwrite_u64(Wi_1,64)) return 0;
  if (!matwrite_u64(Wi_2,64)) return 0;
  if (!matwrite_u64(T_1,64)) return 0;
  if (!matwrite_u64(tmp,64)) return 0;
  if (!matwrite_u64(U_1,64)) return 0;
  if (!matwrite_u64(tmp2,64)) return 0;
  if (sizeof(int) != sizeof(s32)) {
    s32 i, sbuf[64];
    for (i = 0; i < 64; i++)
      sbuf[i] = Si[i];
    if (!matwrite_s32(sbuf,64)) return 0;
    for (i = 0; i < 64; i++)
      sbuf[i] = Si_1[i];
    if (!matwrite_s32(sbuf,64)) return 0;
  }
  else {
    if (!matwrite_s32(Si,64)) return 0;
    if (!matwrite_s32(Si_1,64)) return 0;
  }
  if (!matwrite_u64(X,n)) return 0;
  if (!matwrite_u64(Y,n)) return 0;
  if (!matwrite_u64(Vi,n)) return 0;
  if (!matwrite_u64(Vi_1,n)) return 0;
  if (!matwrite_u64(Vi_2,n)) return 0;
  fclose(file);
  return 1;
}
