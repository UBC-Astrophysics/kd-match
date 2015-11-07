/*
This file is part of ``kd-match'', a suite of programs for matching
stellar catalogues with coordinate systems that differ through affine
transformations (rotations, translations, shearing and scaling).
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

int
main(int argc, char *argv[])
{
  unsigned int n1, j;
  unsigned int cols[]={1,2}, ncolumns=2;
  char **argptr, *filename=NULL;
  double pos[2], transform[6], dumx, dumy, ra_c,dec_c;
#define MAXCOLUMNS 100
  char *fs, *inputstring, **ap, *argv2[MAXCOLUMNS];
  FILE *in;
  char buffer[1024];
  int dotransform=0, doimage=0;

  if (argc<2) {
    printf("Format:\n\n   transform file1\n\n\
   where the options can appear anywhere in any order:\n\n\
   -x column                 column to read x-coordinate (or RA in degrees) from file - default %d\n\
   -y column                 column to read y-coordinate (or Dec in degrees) from file - default %d\n\
   -fs FS                    field separator - default space/TAB\n\
   -c  RA_Centre Dec_Centre  set field center and use x and y as RA and Dec\n\
   -t  params                six parameter transformation to apply (from triangle_kd)\n\
   -                         read from standard input\n\n\
   Only the first files listed will be read.  The final listed parameter\n\
   stands.\n\n\
   transform transforms one set of coordinates to another.  If the field\n\
   center is given, it will return the distance north and east of center\n\
   in arcseconds as x and y.  The transformation will act on these coordinates.\n",cols[0],cols[1]);
    return -1;
  }

  fs= strdup(" \t");

  for (argptr=argv+1;argptr<argv+argc;argptr++) {
    if (strstr(*argptr,"-x")) {
      if (++argptr<argv+argc) {
	cols[0]=atoi(*argptr);
      }
    } else if (strstr(*argptr,"-y")) {
      if (++argptr<argv+argc) {
	cols[1]=atoi(*argptr);
      }
    } else if (strstr(*argptr,"-fs")) {
      if (++argptr<argv+argc) {
	free ( (void *) fs);
	fs=strdup(*argptr);
      }
    } else if (strstr(*argptr,"-t")) {
      dotransform=1;
      for (j=0;j<6 && ++argptr<argv+argc;j++) {
	transform[j]=atof(*argptr);
      }
    } else if (strstr(*argptr,"-c")) {
      doimage=1;
      if (++argptr<argv+argc) {
	ra_c=atof(*argptr)*M_PI/180.0;
      }
      if (++argptr<argv+argc) {
	dec_c=atof(*argptr)*M_PI/180.0;
      }
    } else {
      if (filename==NULL) 
	filename=*argptr;
    }
  }
  /* now read in catalogue 1 */
  if (strcmp(filename,"-")) {
    if ((in=fopen(filename,"r"))==NULL) {
      printf("Unable to open %s  __FILE__:__LINENO__\n",filename);
      return 0;
    }
  } else {
    in=stdin;
  }

  while (fgets(buffer,1023,in)) {
    if (buffer[0]=='#') {
      fputs(buffer,stdout);
    } else {
      inputstring=strdup(buffer);
      /* break line into up to MAXCOLUMNS columns */
      for (ap = argv2; (*ap = strsep(&inputstring, fs)) != NULL;)
	if (**ap != '\0')
	  if (++ap >= &argv2[MAXCOLUMNS])
	    break;
      /* were there any tokens? */
      if (ap>argv2) {
	/* assign columns to the data arrays; missing values given nan */
	for (j=0;j<ncolumns;j++) {
	  pos[j]=(argv2+cols[j]<=ap ? atof(argv2[cols[j]-1]) : 0.0/0.0);
	}
	if (doimage) {
	  pos[0]*=M_PI/180.0; pos[1]*=M_PI/180.0;
	  dumy=acos(cos(pos[1])*cos(dec_c)*cos(pos[0]-ra_c)+sin(pos[1])*sin(dec_c));
	  dumx=atan2(sin(pos[1])-cos(dumy)*sin(dec_c),cos(pos[1])*sin(pos[0]-ra_c)*cos(dec_c));
	  dumy*=180.0*3600.0/M_PI;
	  pos[0]=dumy*cos(dumx);
	  pos[1]=dumy*sin(dumx);
	}
	if (dotransform) {
	  dumx=pos[0]*transform[0]+pos[1]*transform[1]+transform[2];
	  pos[1]=pos[0]*transform[3]+pos[1]*transform[4]+transform[5];
	  pos[0]=dumx;
	}
	if (!isnan(pos[0]) && !isnan(pos[1]))
	  printf("%8.4f %8.4f %s",pos[0],pos[1],buffer);
      }
      free((void *) inputstring);
    }
  }
  if (in!=stdin) fclose(in);
  free ( (void *) fs);


}
