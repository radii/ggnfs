/*********************************************************************
 * matsave.c                                                         *
 * Save matsolve state at regular intervals and when interrupted.    *
 * Copyright 2006, G W Reynolds                                      *
 * TODO:                                                             *
 * If there is an error in the save file, resume from the backup.    *
 * Make stronger checks that the save file matches the loaded matrix.*
 *********************************************************************/

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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
 *   02111-1307  USA
*/

#include <signal.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include "ggnfs.h"
#include "if.h"

/* /etc/magic entry:

   0       ulong   0x4040513B      GGNFS matsolve save file
   >4      ulong   x               version %ld,
   >8      long    x               %ld columns,
   >12     ulong   x               iteration %lu.
*/
#define MATSAVE_MAGIC 0x4040513B
#define MATSAVE_VERSION 1
#define MATSAVE_FILE_NAME "matsave"
#define MATSAVE_BACKUP_NAME "matsave.bak"

/* Exported variables */
volatile int save_flag = 0;
volatile int quit_flag = 0;
int matsave_interval = 0;

#ifdef __MINGW32__
 #ifndef SIGALRM
  #define SIGALRM 14
 #endif

 #include <windows.h>
 void CALLBACK raise_SIGALRM(HWND hwnd, UINT uMsg, UINT_PTR idEvent, DWORD dwTime);
 #define alarm(delay) SetTimer((HWND)NULL, (UINT)NULL, delay, (TIMERPROC) raise_SIGALRM)

 void CALLBACK raise_SIGALRM(HWND hwnd, UINT uMsg, UINT_PTR idEvent, DWORD dwTime) {
    raise(SIGALRM);
 }
#endif

/* Local variables */
static u32 matsave_magic = MATSAVE_MAGIC;
static u32 matsave_version = MATSAVE_VERSION;
static FILE *save_file = NULL;
static int matsave_error = 0;

static void handle_signal(int signum)
{
  switch (signum) {
    case SIGALRM:
      save_flag = 1;
      alarm(matsave_interval); /* reset the alarm */
      break;
    case SIGTERM:
    case SIGINT:
      save_flag = 1;
      quit_flag = 1;
      break;
  }
}

/* Set up signal handlers. */
static void init_matsave(void)
{
  if (matsave_interval > 0) {
    signal(SIGALRM, handle_signal);
    alarm(matsave_interval);
  }
  signal(SIGINT, handle_signal);
  signal(SIGTERM, handle_signal);
}

static int errorclose0(const char *m)
{
  matsave_error = 1;
  fprintf(stderr, m);
  fclose(save_file);
  return 0;
}

static int matread_u32(u32 *d, size_t s)
{
  if (read_u32(save_file, d, s) != s)
    return errorclose0("Error while reading from save file.\n");
  return 1;
}

static int matread_s32(s32 *d, size_t s)
{
  if (read_i32(save_file, d, s) != s)
    return errorclose0("Error while reading from save file.\n");
  return 1;
}

static int matread_u64(u64 *d, size_t s)
{
  if (read_u64(save_file, d, s) != s)
    return errorclose0("Error while reading from save file.\n");
  return 1;
}

/* Attempt to resume from save file. Return the current iteration if
   successful or 0 if not. */
u32 matresume(s32 n, u64 *Wi, u64 *Wi_1, u64 *Wi_2, u64 *T, u64 *T_1,
              u64 *tmp, u64 *U_1, u64 *tmp2, int *Si, int *Si_1,
              u64 *X, u64 *Y, u64 *Vi, u64 *V0, u64 *Vi_1, u64 *Vi_2,
              u64 *tmp_n,u64 *tmp2_n)
{
  u32 magic, version, iterations;
  s32 columns, i, sbuf[64];

  init_matsave();
  if ((save_file = fopen(MATSAVE_FILE_NAME,"rb")) == NULL) {
    if (errno != ENOENT)
      return errorclose0("Could not open save file for reading.\n");
    else
      return 0; /* No save file present */
  }
  if (!matread_u32(&magic,1)) return 0;
  if (magic != matsave_magic)
    return errorclose0("Wrong magic number in save file.\n");
  if (!matread_u32(&version,1)) return 0;
  if (version > matsave_version)
    return errorclose0("Unknown save file version.\n");
  if (!matread_s32(&columns,1)) return 0;
  if (n != columns)
    return errorclose0("Save file does not match loaded matrix.\n");
  if (!matread_u32(&iterations,1)) return 0;
  if (!matread_u64(Wi,64)) return 0;
  if (!matread_u64(Wi_1,64)) return 0;
  if (!matread_u64(Wi_2,64)) return 0;
  if (!matread_u64(T,64)) return 0;
  if (!matread_u64(T_1,64)) return 0;
  if (!matread_u64(tmp,64)) return 0;
  if (!matread_u64(U_1,64)) return 0;
  if (!matread_u64(tmp2,64)) return 0;
  /* Si and Si_1 are arrays of int, but are stored as arrays of s32 so
     that the save file is architecture independent. */
  if (!matread_s32(sbuf,64)) return 0;
  for (i = 0; i < 64; i++)
    Si[i] = sbuf[i];
  if (!matread_s32(sbuf,64)) return 0;
  for (i = 0; i < 64; i++)
    Si_1[i] = sbuf[i];
  if (!matread_u64(X,n)) return 0;
  if (!matread_u64(Y,n)) return 0;
  if (!matread_u64(Vi,n)) return 0;
  if (!matread_u64(V0,n)) return 0;
  if (!matread_u64(Vi_1,n)) return 0;
  if (!matread_u64(Vi_2,n)) return 0;
  if (!matread_u64(tmp_n,n)) return 0;
  if (!matread_u64(tmp2_n,n)) return 0;
  fclose(save_file);
  fprintf(stdout, "Resuming from save file at iteration %" PRIu32 ".\n",
          iterations);

  return iterations;
}

static int matwrite_u32(u32 *d, size_t s)
{
  if (write_u32(save_file, d, s) != s)
    return errorclose0("Error while writing to save file.\n");
  return 1;
}

static int matwrite_s32(s32 *d, size_t s)
{
  if (write_i32(save_file, d, s) != s)
    return errorclose0("Error while writing to save file.\n");
  return 1;
}

static int matwrite_u64(u64 *d, size_t s)
{
  if (write_u64(save_file, d, s) != s)
    return errorclose0("Error while writing to save file.\n");
  return 1;
}

/* Attempt to create a save file. Return 1 if successful, 0 if not. */
int matsave(u32 iterations, s32 n, u64 *Wi, u64 *Wi_1, u64 *Wi_2, u64 *T,
            u64 *T_1, u64 *tmp, u64 *U_1, u64 *tmp2, int *Si, int *Si_1,
            u64 *X, u64 *Y, u64 *Vi, u64 *V0, u64 *Vi_1, u64 *Vi_2,
            u64 *tmp_n, u64 *tmp2_n)
{
  s32 i, sbuf[64];

  fprintf(stdout, "\nCreating save file at iteration %" PRIu32 ".\n",
          iterations);
  /* Don't overwrite backup if resume or an earlier save had problems */
  if (!matsave_error && rename(MATSAVE_FILE_NAME,MATSAVE_BACKUP_NAME) == -1
      && errno != ENOENT)
    fprintf(stderr, "Could not create backup save file.\n");
  if ((save_file = fopen(MATSAVE_FILE_NAME,"wb")) == NULL)
    return errorclose0("Could not open save file for writing.\n");
  matsave_error = 0;
  if (!matwrite_u32(&matsave_magic,1)) return 0;
  if (!matwrite_u32(&matsave_version,1)) return 0;
  if (!matwrite_s32(&n,1)) return 0;
  if (!matwrite_u32(&iterations,1)) return 0;
  if (!matwrite_u64(Wi,64)) return 0;
  if (!matwrite_u64(Wi_1,64)) return 0;
  if (!matwrite_u64(Wi_2,64)) return 0;
  if (!matwrite_u64(T,64)) return 0;
  if (!matwrite_u64(T_1,64)) return 0;
  if (!matwrite_u64(tmp,64)) return 0;
  if (!matwrite_u64(U_1,64)) return 0;
  if (!matwrite_u64(tmp2,64)) return 0;
  for (i = 0; i < 64; i++)
    sbuf[i] = Si[i];
  if (!matwrite_s32(sbuf,64)) return 0;
  for (i = 0; i < 64; i++)
    sbuf[i] = Si_1[i];
  if (!matwrite_s32(sbuf,64)) return 0;
  if (!matwrite_u64(X,n)) return 0;
  if (!matwrite_u64(Y,n)) return 0;
  if (!matwrite_u64(Vi,n)) return 0;
  if (!matwrite_u64(V0,n)) return 0;
  if (!matwrite_u64(Vi_1,n)) return 0;
  if (!matwrite_u64(Vi_2,n)) return 0;
  if (!matwrite_u64(tmp_n,n)) return 0;
  if (!matwrite_u64(tmp2_n,n)) return 0;
  fclose(save_file);

  return 1;
}
