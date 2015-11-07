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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "kdtree.h"

int verbose=0;

#ifndef __AVAILABILITY__
void
__sincospi(double ang, double *sinval, double *cosval) {
  ang*=M_PI;
  *sinval=sin(ang);
  *cosval=cos(ang);
}
#endif

int
main(int argc, char *argv[]) {
  float dist;
  FILE *in;
  char buffer[1024];
  struct kdtree *kd;
  struct kdres *res;
  float pos[3], irpos[3];
#define MAXCOLUMNS 100
  char *optline, *inputstring, **ap, *argv2[MAXCOLUMNS];
  int cols1[2]={1,2}, cols2[2]={1,2}, ncolumns=2, j;
  int dotransform1=0, dotransform2=0, dounique=0, donearest=1, dosphere=0, loadon=1;
  char *fs1, *fs2, *filename1=NULL, *filename2=NULL;
  double transform1[6], transform2[6], dumx, distance=-10;

  fs1 = strdup(" \t");
  fs2 = strdup(" \t");

  if (argc<3) {
    printf("Format:\n\n   match_kd file1 file2 [options]\n\n\
Find the closest star in catalogue 2 for each star in catalogue 1.\n\
Print the line from catalogue 1 cat2-coords distance the line in catalogue 2.\n\
The options can appear anywhere in any order:\n\n\
   -x1 column    column to read x-coordinate from file 1 - default %d\n\
   -y1 column    column to read y-coordinate from file 1 - default %d\n\
   -x2 column    column to read x-coordinate from file 2 - default %d\n\
   -y2 column    column to read y-coordinate from file 2 - default %d\n\
   -t  params    six parameter transformation from 1 to 2 (from triangle_kd)\n\
   -t2 params    six parameter transformation from 2 to 1 (from triangle_kd)\n\
   -d  distance  find all objects in catalogue 2 within the given distance;\n\
                 listed after the closest one (may or may not include the\n\
                 closest object)\n\
   -s            do not output the nearest object\n\
   -n            find all objects in catalogue 2 that are outside the given distance\n\
   -fs  FS       field separator - default space/TAB\n\
   -fs1 FS       field separator for file 1\n\
   -fs2 FS       field separator for file 2\n\
   -eq           coordinates are RA/Dec or l/b on a sphere in degrees\n\
                 (distance here is the areal distance)\n\
   -             read from standard input\n\n\
   Only the first two files listed will be read.  The final listed parameter\n\
   stands.\n\
",cols1[0],cols1[1],cols2[0],cols2[1]);
    return -1;
  }

  printf("#");
  for (ap=argv;ap<argv+argc;ap++) {
    printf(" %s",*ap);
  }
  printf("\n");

  for (ap=argv+1;ap<argv+argc;ap++) {
    if (strstr(*ap,"-x1")) {
      if (++ap<argv+argc) {
	cols1[0]=atoi(*ap);
      }
    } else if (strstr(*ap,"-y1")) {
      if (++ap<argv+argc) {
	cols1[1]=atoi(*ap);
      }
    } else if (strstr(*ap,"-x2")) {
      if (++ap<argv+argc) {
	cols2[0]=atoi(*ap);
      } 
    } else if (strstr(*ap,"-y2")) {
      if (++ap<argv+argc) {
	cols2[1]=atoi(*ap);
      }
    } else if (strstr(*ap,"-fs1")) {
      if (++ap<argv+argc) {
	free ( (void *) fs1);
	fs1=strdup(*ap);
      }
    } else if (strstr(*ap,"-fs2")) {
      if (++ap<argv+argc) {
	free ( (void *) fs2);
	fs2=strdup(*ap);
      }
    } else if (strstr(*ap,"-fs")) {
      if (++ap<argv+argc) {
	free ( (void *) fs1);
	fs1=strdup(*ap);
	free ( (void *) fs2);
	fs2=strdup(*ap);
      }
    } else if (strstr(*ap,"-d")) {
      distance=atof(*(++ap));
    } else if (strstr(*ap,"-s")) {
      donearest=0;
    } else if (strstr(*ap,"-n")) {
      dounique=1;
    } else if (strstr(*ap,"-eq")) {
      dosphere=1;
    } else if (strstr(*ap,"-t2")) {
      dotransform2=1;
      for (j=0;j<6 && ++ap<argv+argc;j++) {
	transform2[j]=atof(*ap);
      }
    } else if (strstr(*ap,"-t")) {
      dotransform1=1;
      for (j=0;j<6 && ++ap<argv+argc;j++) {
	transform1[j]=atof(*ap);
      }
    } else if (strstr(*ap,"-v")) {
      verbose++;
    } else {
      if (filename1==NULL) 
	filename1=*ap;
      else if (filename2==NULL) 
	filename2=*ap;
    }
  }

  if (verbose) {
    printf("#");
    for (j=0;j<argc;j++) {
      printf(" %s",argv[j]);
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
    if (dotransform1) {
      printf("# transform1=");
      for (j=0;j<6;j++) {
	printf("%g ",transform1[j]);
      }
      printf("\n");
    }
    if (dotransform2) {
      printf("# transform2=");
      for (j=0;j<6;j++) {
	printf("%g ",transform2[j]);
      }
      printf("\n");
    }
  }

  /* Find the closest star in catalogue 2 for each star in catalogue 1 */
  /* Print the line from catalogue 2 - distance - the line in catalogue 1 */

  /* create the kd-tree for the star positions */
  if (dosphere) {
    kd = kd_create(2);
  } else {
    kd = kd_create(3);
  }
  /* designate a function to deallocate the data */
  kd_data_destructor(kd,free);

  /* actually read in catalogue 2 first */
  if (strcmp(filename2,"-")) {
    if ((in=fopen(filename2,"r"))==NULL) {
      printf("Unable to open %s  __FILE__:__LINENO__\n",filename2);
      return 0;
    }
  } else {
    in=stdin;
  }
  loadon=1;
  while (fgets(buffer,1023,in)) {
    if (buffer[0]=='*') loadon=1-loadon;
    if (buffer[0]=='#' || buffer[0]=='*' || !loadon) {
      fputs(buffer,stdout);
    } else {
      optline=strdup(buffer);
      inputstring=buffer;
      /* break line into up to MAXCOLUMNS columns */
      for (ap = argv2; (*ap = strsep(&inputstring, fs2)) != NULL;)
	if (**ap != '\0')
	  if (++ap >= &argv2[MAXCOLUMNS])
	    break;
      /* were there any tokens? */
      if (ap>argv2) {
	/* assign columns to the data arrays; missing values given nan */
	for (j=0;j<ncolumns;j++) {
	  pos[j]=(argv2+cols2[j]<=ap ? atof(argv2[cols2[j]-1]) : 0.0/0.0);
	}
      }
      if (dotransform2) {
	  /* printf("\n%g %g ",irpos[0],irpos[1]);*/
	  dumx=pos[0]*transform2[0]+pos[1]*transform2[1]+transform2[2];
	  pos[1]=pos[0]*transform2[3]+pos[1]*transform2[4]+transform2[5];
	  pos[0]=dumx;
	  /* printf("%g %g\n",irpos[0],irpos[1]);*/
      }
      if (dosphere) {
	double x,y,z,rad;
	__sincospi(pos[0]/180.0,&y,&x);
	__sincospi(pos[1]/180.0,&z,&rad);
	pos[0]=x*rad;
	pos[1]=y*rad;
	pos[2]=z;
      }
      if (kd_insertf(kd, pos, (void *) optline)) {
	printf("Unable to insert point into the tree at  __FILE__:__LINENO__\n");
	return -1;
      }
    }
    /*    printf("%ld %g %g\n",xptr-xopt,*(xptr-1),*(yptr-1));  */
  }
  if (in!=stdin) fclose(in);

  /* now read in catalogue 1 */
  if (strcmp(filename1,"-")) {
    if ((in=fopen(filename1,"r"))==NULL) {
      printf("Unable to open %s  __FILE__:__LINENO__\n",filename1);
      return 0;
    }
  } else {
    in=stdin;
  }
  loadon=1;
  while (fgets(buffer,1023,in)) {
    if (buffer[0]=='*') loadon=1-loadon;
    if (buffer[0]=='#' || buffer[0]=='*' || !loadon) {
      fputs(buffer,stdout);
    } else {
      buffer[strlen(buffer)-1]=0;
      if (!dounique && donearest) {
	printf("%s",buffer);
      }
      inputstring=strdup(buffer);
      /* break line into up to MAXCOLUMNS columns */
      for (ap = argv2; (*ap = strsep(&inputstring, fs1)) != NULL;)
	if (**ap != '\0')
	  if (++ap >= &argv2[MAXCOLUMNS])
	    break;
      /* were there any tokens? */
      if (ap>argv2) {
	/* assign columns to the data arrays; missing values given nan */
	for (j=0;j<ncolumns;j++) {
	  irpos[j]=(argv2+cols1[j]<=ap ? atof(argv2[cols1[j]-1]) : 0.0/0.0);
	}
	if (dotransform1) {
	  /* printf("\n%g %g ",irpos[0],irpos[1]);*/
	  dumx=irpos[0]*transform1[0]+irpos[1]*transform1[1]+transform1[2];
	  irpos[1]=irpos[0]*transform1[3]+irpos[1]*transform1[4]+transform1[5];
	  irpos[0]=dumx;
	  /* printf("%g %g\n",irpos[0],irpos[1]);*/
	}
	if (dosphere) {
	  double x,y,z,rad;
	  __sincospi(irpos[0]/180.0,&y,&x);
	  __sincospi(irpos[1]/180.0,&z,&rad);
	  irpos[0]=x*rad;
	  irpos[1]=y*rad;
	  irpos[2]=z;
	}
	if (!dounique && donearest) {
	  res=kd_nearestf(kd,irpos);
	  if (kd_res_size(res)>0) {
	    optline = (char *) kd_res_itemf( res, pos );
	    dist = hypot(pos[0]-irpos[0],pos[1]-irpos[1]);
	    if (dotransform1) {
	      printf(" %8.4f %8.4f",irpos[0],irpos[1]);
	    }
	    if (dotransform2) {
	      printf(" %8.4f %8.4f",pos[0],pos[1]);
	    }
	    printf(" %12.4e %s",dist,optline);
	  }
	}
	if (distance>0) {
	  res=kd_nearest_rangef(kd,irpos,distance);
	  if (kd_res_size(res)>0) {
	    while( !kd_res_end( res ) ) {
	      optline = (char *) kd_res_itemf( res, pos );
	      if (dounique) {
		optline[0]=0;
	      } else {
		if (dosphere) {
		  dist = hypot(hypot(pos[0]-irpos[0],pos[1]-irpos[1]),pos[2]-irpos[2]);
		} else {
		  dist = hypot(pos[0]-irpos[0],pos[1]-irpos[1]);
		}
		if (dotransform1) {
		  printf("%s %8.4f %8.4f %8.4f %s",buffer,irpos[0],irpos[1],dist,optline);
		} else {
		  printf("%s %8.4f %s",buffer,dist,optline);
		}
	      }
	      /* go to the next entry */
	      kd_res_next( res );
	    }
	  }
	}
      }
    }
    free((void *) inputstring);
  }
  if (in!=stdin) fclose(in);

  if (dounique) {
    /* print out all of the stars in catalogue 2 */
    /* that were outside the search radius */ 
    irpos[0]=0; irpos[1]=0;
    res=kd_nearest_rangef(kd,irpos,1e100);
    if (kd_res_size(res)>0) {
      while( !kd_res_end( res ) ) {
	optline = (char *) kd_res_itemf( res, pos );
	if (optline[0]) {
	  printf("%s",optline);
	}
	/* go to the next entry */
	kd_res_next( res );
      }
    }
  }

  kd_free(kd);
  free ( (void *) fs1);
  free ( (void *) fs2);

  return(0);
}
