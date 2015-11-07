/*
This file is part of ``kd-match'', a suite of programs for matching stellar catalogues
with coordinate systems that differ through affine transformations (rotations, 
translations, shearing and scaling). 
Copyright (C) 2013 Jeremy Heyl <heyl@phas.ubc.ca>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "loadfile.h"
#include "calctransform.h"

int
main(int argc, char *argv[])
{
  unsigned int n1, n2, np;
  unsigned int cols1[]={1,2}, cols2[]={1,2};
  char **argptr, *filename1=NULL, *filename2=NULL;
  double *d1ptr[2], *d2ptr[2], paramx[3], paramy[3], x, y;
  char *fs1, *fs2;
  int verbose, i, dotransform=0;

  fs1 = strdup(" \t");
  fs2 = strdup(" \t");

  if (argc<3) {
    printf("Format:\n\n   calctrans file1 file2 [options]\n\n\
   where the options can appear anywhere in any order:\n\n\
   -x1 column          column to read x-coordinate from file 1 - default %d\n\
   -y1 column          column to read y-coordinate from file 1 - default %d\n\
   -x2 column          column to read x-coordinate from file 2 - default %d\n\
   -y2 column          column to read y-coordinate from file 2 - default %d\n\
   -x  x-coord         x-coord to transform\n\
   -y  y-coord         x-coord to transform\n\
   -fs  FS             field separator - default space/TAB\n\
   -fs1 FS             field separator for file 1\n\
   -fs2 FS             field separator for file 2\n\
   -v                  be verbose (more -v more verbose)\n\
   -                   read from standard input\n\n\
   Only the first two files listed will be read.  The final listed parameter\n\
   stands.\n\n\
   calctrans determines the best fitting transformation between the two sets\n\
   of coordinates.  It does not fit asterisms, but rather it assumes that the\n\
   objects in the first list corresponds to those in the second list.\n\n\
",cols1[0],cols1[1],cols2[0],cols2[1]);
    return -1;
  }

  for (argptr=argv+1;argptr<argv+argc;argptr++) {
    if (strstr(*argptr,"-x1")) {
      if (++argptr<argv+argc) {
	cols1[0]=atoi(*argptr);
      }
    } else if (strstr(*argptr,"-y1")) {
      if (++argptr<argv+argc) {
	cols1[1]=atoi(*argptr);
      }
    } else if (strstr(*argptr,"-x2")) {
      if (++argptr<argv+argc) {
	cols2[0]=atoi(*argptr);
      } 
    } else if (strstr(*argptr,"-y2")) {
      if (++argptr<argv+argc) {
	cols2[1]=atoi(*argptr);
      }
    } else if (strstr(*argptr,"-x")) {
      if (++argptr<argv+argc) {
	dotransform=1;
	x=atof(*argptr);
      }
    } else if (strstr(*argptr,"-y")) {
      if (++argptr<argv+argc) {
	dotransform=1;
	y=atof(*argptr);
      }
    } else if (strstr(*argptr,"-fs1")) {
      if (++argptr<argv+argc) {
	free ( (void *) fs1);
	fs1=strdup(*argptr);
      }
    } else if (strstr(*argptr,"-fs2")) {
      if (++argptr<argv+argc) {
	free ( (void *) fs2);
	fs2=strdup(*argptr);
      }
    } else if (strstr(*argptr,"-fs")) {
      if (++argptr<argv+argc) {
	free ( (void *) fs1);
	fs1=strdup(*argptr);
	free ( (void *) fs2);
	fs2=strdup(*argptr);
      }
    } else if (strstr(*argptr,"-v")) {
      verbose++;
    } else {
      if (filename1==NULL) 
	filename1=*argptr;
      else if (filename2==NULL) 
	filename2=*argptr;
    }
  }

  if (verbose) {
    printf("#");
    for (i=0;i<argc;i++) {
      printf(" %s",argv[i]);
    }
    printf("\n");
  }

  if (verbose>1) {
    printf("# xcol1= %d\n",cols1[0]);
    printf("# ycol1= %d\n",cols1[1]);
    printf("# xcol2= %d\n",cols2[0]);
    printf("# ycol2= %d\n",cols2[1]);
    printf("# fs1= %s\n",fs1);
    printf("# fs2= %s\n",fs2);
    printf("# %s %s\n",filename1,filename2);
  }

  if (strcmp(filename1,"-")) {
    n1=loadfile_fs(filename1,2,cols1,d1ptr,fs1);
  } else {
    n1=loadfile_fileptr_fs(stdin,2,cols1,d1ptr,fs1);
  }

  if (filename2==NULL)  { filename2=filename1; }

  if (strcmp(filename2,"-")) {
    n2=loadfile_fs(filename2,2,cols2,d2ptr,fs2); 
  } else {
    n2=loadfile_fileptr_fs(stdin,2,cols2,d2ptr,fs2); 
  }

  np = (n1<n2 ? n1 : n2);
  calctransform_noindex(d1ptr[0],d1ptr[1],d2ptr[0],np,paramx);
  printf("Transforms from 1->2\n");
  printf("x2 = %g x1 + %g y1 + %g\n",paramx[0],paramx[1],paramx[2]);

  if (dotransform) { printf("x2 = %g\n",paramx[0]*x+paramx[1]*y+paramx[2]); }
  calctransform_noindex(d1ptr[0],d1ptr[1],d2ptr[1],np,paramy);
  printf("y2 = %g x1 + %g y1 + %g\n",paramy[0],paramy[1],paramy[2]);
  if (dotransform) { printf("y2 = %g\n",paramy[0]*x+paramy[1]*y+paramy[2]); }
  printf("-t %g %g %g %g %g %g\n",paramx[0],paramx[1],paramx[2],paramy[0],paramy[1],paramy[2]);

  printf("Transforms from 2->1\n");
  calctransform_noindex(d2ptr[0],d2ptr[1],d1ptr[0],np,paramx);
  printf("x1 = %g x2 + %g y2 + %g\n",paramx[0],paramx[1],paramx[2]);
  if (dotransform) { printf("x1 = %g\n",paramx[0]*x+paramx[1]*y+paramx[2]); }
  calctransform_noindex(d2ptr[0],d2ptr[1],d1ptr[1],np,paramy);
  printf("y1 = %g x2 + %g y2 + %g\n",paramy[0],paramy[1],paramy[2]);
  if (dotransform) { printf("y1 = %g\n",paramy[0]*x+paramy[1]*y+paramy[2]); }
  printf("-t %g %g %g %g %g %g\n",paramx[0],paramx[1],paramx[2],paramy[0],paramy[1],paramy[2]);

  free((void *) d2ptr[0]);
  free((void *) d2ptr[1]);
  free((void *) d1ptr[0]);
  free((void *) d1ptr[1]);

  return 0;

}
