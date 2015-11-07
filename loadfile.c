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

#define ALLOCBLOCK 512
#define MAXCOLUMNS 100
#define LOADFILE_FS " \t"

unsigned int
loadfile_fileptr_fs(FILE *in, unsigned int ncolumns, unsigned int columns[], double *data[], char *fs) {
  unsigned int i, j, ialloc;
  char buffer[1024];
  char **ap, *argv[MAXCOLUMNS], *inputstring;
  int loadon;

  for (j=0;j<ncolumns;j++) {
    if ((data[j]=(double *) malloc(sizeof(double)*512))==NULL) {
      printf("Unable to allocate data[%d] in __FILE__:__LINENO__\n",j);
      return 0;
    }
  }
  ialloc=512;

  i=0;
  loadon=1;
  while (fgets(buffer,1023,in)) {
    if (buffer[0]=='*') loadon=1-loadon;
    if (buffer[0]=='#' || buffer[0]=='*' || !loadon) {
      /* skip the line */
    } else {
      inputstring=buffer;
      /* break line into up to MAXCOLUMNS columns */
      for (ap = argv; (*ap = strsep(&inputstring, fs)) != NULL;)
	if (**ap != '\0')
	  if (++ap >= &argv[MAXCOLUMNS])
	    break;
      /* were there any tokens? */
      if (ap>argv) {
	/* assign columns to the data arrays; missing values given nan */
	for (j=0;j<ncolumns;j++) {
	  data[j][i]=(argv+columns[j]<=ap ? atof(argv[columns[j]-1]) : 0.0/0.0);
	}
	i++;
	/* do we need to allocate more memory? */
	if (i==ialloc) {
	  for (j=0;j<ncolumns;j++) {
	    if ((data[j]=(double *) realloc((void *) data[j],sizeof(double)*(i+512)))==NULL) {
	      printf("Unable to reallocate data[%d] in __FILE__:__LINENO__\n",j);
	      return 0;
	    }
	  }
	  ialloc=i+512;
	}
      }
    }
  }
  /* free excess memory */
  for (j=0;j<ncolumns;j++) {
    if ((data[j]=(double *) realloc((void *) data[j],sizeof(double)*i))==NULL) {
      printf("Unable to reallocate data[%d] in __FILE__:__LINENO__\n",j);
      return 0;
    }
  }
  
  /* return the number of lines read */
  return i;
}

unsigned int
loadfile_fs(char *filename, unsigned int ncolumns, unsigned int columns[], double *data[], char *fs) {
  FILE *in;
  unsigned int retval;

  if ((in=fopen(filename,"r"))==NULL) {
    printf("Unable to open %s  __FILE__:__LINENO__\n",filename);
    return 0;
  }

  retval=loadfile_fileptr_fs(in, ncolumns, columns, data, fs);
  fclose(in);
  return retval;
}

unsigned int
loadfile_fileptr(FILE *in, unsigned int ncolumns, unsigned int columns[], double *data[]) {
  return loadfile_fileptr_fs(in,ncolumns,columns,data,LOADFILE_FS);
}

unsigned int
loadfile(char *filename, unsigned int ncolumns, unsigned int columns[], double *data[]) {
  return loadfile_fs(filename,ncolumns,columns,data,LOADFILE_FS);
}
