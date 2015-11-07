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
#include <assert.h>
#include "kdtree.h"
#include "loadfile.h"


unsigned int n1, n2;
double *xp1, *yp1, *xp2, *yp2;
int listswapped=0, verbose=0, noswap=0, matching_pairs;
struct kdtree *kd_good;
double dist_cut=0.2, trans_cut=0.2, x1_factor=1, y1_factor=1;
double bestcoeff[2];
int nbest=0;


void
calctranslate(double xp1[],double yp1[],double xp2[],double yp2[],
	      int index1[],int index2[],int npoint,double param[]) {
  int i;
  param[0]=param[1]=0;
  for (i=0;i<npoint;i++) {
    param[0]+=xp2[index2[i]]-xp1[index1[i]];
    param[1]+=yp2[index2[i]]-yp1[index1[i]];
  }
  param[0]/=npoint;
  param[1]/=npoint;
}


int
addpair(double *param, int *pairdata) {
  struct kdres *res;
  static int pair_added;
  double pos[2], coeff[2];
  int *data, i;

  matching_pairs=0;
  if (pair_added) {
    res=kd_nearest_range(kd_good,param,trans_cut);
    /* if there are some pairs with matching transforms, then tell us about them */
    if (kd_res_size(res)>0) {
      double sx1, sy1, s33;
#if 0
      sx1=sy1=s33=0;
      data=pairdata;
      for (i=0;i<2;i++) {
	if (listswapped) {
	  sx1+=xp1[data[i+2]]-xp2[data[i]];
	  sy1+=yp1[data[i+2]]-yp2[data[i]];
	} else {
	  sx1+=xp2[data[i+2]]-xp1[data[i]];
	  sy1+=yp2[data[i+2]]-yp1[data[i]];
	}
	s33++;
      }
      s33++;
#endif
      sx1=sy1=s33=0;

      matching_pairs=kd_res_size(res);
      if (verbose>=0) {
	printf("x-transform: x2= %g x1 + 0 y1 + %g\n",x1_factor,param[0]);
	printf("y-transform: y2= 0 x1 + %g y1 + %g\n",y1_factor,param[1]);
	printf("Number of matching pairs: %d\n",matching_pairs);
      } 
      while( !kd_res_end( res ) ) {
	/* get the data and position of the current result item */
	data = (int*)kd_res_item( res, pos );
	if (verbose>=0) {
	  printf("Pair with matching translation %g %g {",pos[0],pos[1]);
	  for (i=0;i<2;i++) {
	    printf(" %d",data[i]);
	  }
	  printf("} -> {");
	  for (i=2;i<4;i++) {
	    printf(" %d",data[i]);
	  }
	  printf("(%g - %g)",xp1[data[0]],xp2[data[2]]);
	  printf("}\n");
	}
	for (i=0;i<2;i++) {
	  if (listswapped) {
	    sx1+=xp1[data[i+2]]-xp2[data[i]];
	    sy1+=yp1[data[i+2]]-yp2[data[i]];
	  } else {
	    sx1+=xp2[data[i+2]]-xp1[data[i]];
	    sy1+=yp2[data[i+2]]-yp1[data[i]];
	  }
	  s33++;
	}
	/* go to the next entry */
	kd_res_next( res );
      }
      coeff[0]=sx1/s33;
      coeff[1]=sy1/s33;
      if (verbose>=0) {
	printf("Sum Calc:\nx2 = %g x1 + %g y1 + %g\n",
	       x1_factor,0.0,coeff[0]);
	printf("y2 = %g x1 + %g y1 + %g\n -t %g 0 %g 0 %g %g\n",
	       0.0,y1_factor,coeff[1],x1_factor,coeff[0],y1_factor,coeff[1]);
      }
      if (matching_pairs>nbest) {
	nbest=matching_pairs;
	for (i=0;i<2;i++) {
	  bestcoeff[i]=coeff[i];
	}
      }
    }
    /* free the results structure */
    kd_res_free(res);
  }
  pair_added=1;

  return kd_insert(kd_good, param, (void *) pairdata);
}

void
pairoutput(int ih, int jh, int iah, int jah) {
  int i;
  int index1[2], index2[2], *data;
  double param[2];

  index1[0]=ih; index1[1]=jh;
  index2[0]=iah; index2[1]=jah;
  if (listswapped) { 
    /* output the pair information */
    if (verbose>1) {
      printf("Pair %3d - %3d: %10.4f %10.4f; %3d - %3d: %10.4f %10.4f\n",
	     iah,jah,xp2[iah]-xp2[jah],yp2[iah]-yp2[jah],
	     ih,jh,xp1[ih]-xp1[jh],yp1[ih]-yp1[jh]);
    }

    /* calculate forward transformation */
    calctranslate(xp2,yp2,xp1,yp1,index2,index1,2,param);
    data=(int *) malloc(sizeof(int)*4);
    for (i=0;i<2;i++) {
      data[i]=index2[i]; 
      data[i+2]=index1[i];
    }

    /* add it to the tree */
    assert(addpair(param,data) == 0);

    if (verbose) {
      printf("Transforms from 1->2\n");
      printf("x2 = %g x1 + 0 y1 + %g\n",x1_factor,param[0]);
      printf("y2 = 0 x1 + %g y1 + %g\n",y1_factor,param[1]);
    }
    
    /* output reverse transformation */
    if (verbose>1) {
      printf("Transforms from 2->1\n");
      calctranslate(xp1,yp1,xp2,yp2,index1,index2,2,param);
      printf("x2 = 1 x1 + 0 y1 + %g\n",param[0]);
      printf("y2 = 0 x1 + 1 y1 + %g\n",param[1]);
    }
  } else {
    /* output the pair information */
    if (verbose>1) {
      printf("Pair %3d - %3d - %10.4f %10.4f - %10.4f %10.4f; %3d - %3d %10.4f %10.4f\n",
	     ih,jh,xp2[ih],xp2[jh],xp2[ih]-xp2[jh],yp2[ih]-yp2[jh],
	     iah,jah,xp1[iah]-xp1[jah],yp1[iah]-yp1[jah]);
    }

    /* calculate forward transformation */
    calctranslate(xp1,yp1,xp2,yp2,index1,index2,2,param);
    data=(int *) malloc(sizeof(int)*4);
    for (i=0;i<2;i++) {
      data[i]=index1[i]; 
      data[i+2]=index2[i];
    }

    /* add it to the tree */
    assert(addpair(param,data) == 0);

    if (verbose) {
      printf("Transforms from 1->2\n");
      printf("x2 = %g x1 + 0 y1 + %g\n",x1_factor,param[0]);
      printf("y2 = 0 x1 + %g y1 + %g\n",y1_factor,param[1]);
    }
    
    /* output reverse transformation */
    if (verbose>1) {
      printf("Transforms from 2->1\n");
      calctranslate(xp2,yp2,xp1,yp1,index2,index1,2,param);
      printf("x2 = 1 x1 + 0 y1 + %g\n",param[0]);
      printf("y2 = 0 x1 + 1 y1 + %g\n",param[1]);
    }
  }
}



int
main(int argc, char *argv[])
{
  int i, j, k, max_matches=20;
  int *data, datahold[2];
  double la[2], diff, pos[2],  atof();
  struct kdtree *kd;
  struct kdres *res;
  unsigned int cols1[]={1,2}, cols2[]={1,2};
  char **argptr, *filename1=NULL, *filename2=NULL;
  double *dptr[2];
  char *fs1, *fs2;

  fs1 = strdup(" \t");
  fs2 = strdup(" \t");

  if (argc<3) {
    printf("Format:\n\n   pair_kd file1 file2 [options]\n\n\
   where the options can appear anywhere in any order:\n\n\
   -d distance_cutoff  how small of a distance to call a match - default %g\n\
   -t transform_cutoff how small of a distance to call a match - default %g\n\
   -xf x1_factor       factor to scale the x1 coordinate       - default %g\n\
   -yf y1_factor       factor to scale the y1 coordinate       - default %g\n\
   -m max_matches      number of matching transforms to quit   - default %d\n\
   -x1 column          column to read x-coordinate from file 1 - default %d\n\
   -y1 column          column to read y-coordinate from file 1 - default %d\n\
   -x2 column          column to read x-coordinate from file 2 - default %d\n\
   -y2 column          column to read y-coordinate from file 2 - default %d\n\
   -fs  FS             field separator - default space/TAB\n\
   -fs1 FS             field separator for file 1\n\
   -fs2 FS             field separator for file 2\n\
   -ns                 Do not swap lists, even if the first is larger\n\
   -v                  be more verbose (more -v more verbose)\n\
   -q                  be less verbose (more -q less verbose)\n\
   -                   read from standard input\n\n\
   Only the first two files listed will be read.  The final listed parameter\n\
   stands.\n\n\
   pair_kd will try to match pairs in the two catalogues.  It is useful\n\
   when one expects there to be only translation between the\n\
   two catalogues but no rotation or shearing.  If you know the scaling between the\n\
   two coordinate systems, use the -xf and -yf parameters. If you do expect\n\
   rotation and scaling, use triangle_kd. If you do expect shearing, use quad_kd.\n",
	   dist_cut,trans_cut,x1_factor, y1_factor, max_matches,cols1[0],cols1[1],cols2[0],cols2[1]);
    return -1;
  }

  /* parse command line */
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
    } else if (strstr(*argptr,"-d")) {
      if (++argptr<argv+argc) {
	dist_cut=atof(*argptr);
      }
    } else if (strstr(*argptr,"-t")) {
      if (++argptr<argv+argc) {
	trans_cut=atof(*argptr);
      }
    } else if (strstr(*argptr,"-xf")) {
      if (++argptr<argv+argc) {
        x1_factor=atof(*argptr);
      }
    } else if (strstr(*argptr,"-yf")) {
      if (++argptr<argv+argc) {
        y1_factor=atof(*argptr);
      }
    } else if (strstr(*argptr,"-m")) {
      if (++argptr<argv+argc) {
	max_matches=atoi(*argptr);
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
    } else if (strstr(*argptr,"-ns")) {
      noswap=1;
    } else if (strstr(*argptr,"-v")) {
      verbose++;
    } else if (strstr(*argptr,"-q")) {
      verbose--;
    } else {
      if (filename1==NULL) 
	filename1=*argptr;
      else if (filename2==NULL) 
	filename2=*argptr;
    }
  }

  if (verbose>0) {
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
    printf("# dist_cut= %g\n",dist_cut);
    printf("# fs1= %s\n",fs1);
    printf("# fs2= %s\n",fs2);
    printf("# noswap= %d\n",noswap);
    printf("# %s %s\n",filename1,filename2);
  }

  if (strcmp(filename1,"-")) {
    n1=loadfile_fs(filename1,2,cols1,dptr,fs1);
  } else {
    n1=loadfile_fileptr_fs(stdin,2,cols1,dptr,fs1);
  }
  xp1=dptr[0]; yp1=dptr[1];

  if (strcmp(filename2,"-")) {
    n2=loadfile_fs(filename2,2,cols2,dptr,fs2); 
  } else {
    n2=loadfile_fileptr_fs(stdin,2,cols2,dptr,fs2); 
  }
  xp2=dptr[0]; yp2=dptr[1];


  if (x1_factor!=1) {
    for (i=0;i<n1;i++) {
      xp1[i]*=x1_factor;
    }
  }

  if (y1_factor!=1) {
    for (i=0;i<n1;i++) {
      yp1[i]*=y1_factor;
    }
  }

  /* output some information at the top */
  if (verbose>0) {
    printf("# List 1: %s has %u points and %lu pairs\n# List 2: %s has %u points and %lu pairs\n",
	   filename1,n1,n1*(n1-1ul)*(n1-2ul)/6,filename2,n2,n2*(n2-1ul)*(n2-2ul)/6);
  }

  /* always build the tree from the smaller list */
  /* it saves memory and time */
  if (n1>n2 && !noswap) {
    /* if the first list is larger than swap the lists */
    double *dumptr;
    dumptr=xp1; xp1=xp2; xp2=dumptr;
    dumptr=yp1; yp1=yp2; yp2=dumptr;
    i=n1; n1=n2; n2=i;
    listswapped=1;
    if (verbose>0) printf("# For internal use the lists were swapped.\n");
  } else {
    listswapped=0;
  }


  /* create the kd-tree for the pairs*/
  kd = kd_create(2);
  /* designate a function to deallocate the data */
  kd_data_destructor(kd,free);

  /* build the tree */
  for (i=0;i<n1-1;i++) {
    for (j=i+1;j<n1;j++) {
	/* calculate difference */
      la[0]=xp1[i]-xp1[j];
      la[1]=yp1[i]-yp1[j];
      /* if x-difference is less than a pixel, skip this pair */
      if (fabs(la[0])<1) break;

      /* allocate an array to hold the points */
      data=(int *) malloc(sizeof(int)*2);
      if (la[0]>0) { la[0]*=-1; la[1]*=-1;
	data[0]=j; data[1]=i; 
      } else {
	data[0]=i; data[1]=j; 
      }
      /* add it to the tree */
      assert(kd_insert(kd, la, (void *) data) == 0);
    }
  }

  /* go through the pairs from the longer list */

  /* create the kd-tree for the matches */
  kd_good = kd_create(2);
  /* designate a function to deallocate the data */
  kd_data_destructor(kd_good,free);

  for (i=0;i<n2-1;i++) {
    for (j=i+1;j<n2;j++) {
	/* calculate side lengths */
        la[0]=xp2[i]-xp2[j];
        la[1]=yp2[i]-yp2[j];

	/* if x-difference is less than a pixel, skip this pair */
	if (fabs(la[0])<1) break;

	if (la[0]>0) { 
	  la[0]*=-1; la[1]*=-1;
	  datahold[0]=j; datahold[1]=i; 
	} else {
	  datahold[0]=i; datahold[1]=j; 
	}
	/* find all the pairs from the first (short) list that are within the dist_cut */
	res=kd_nearest_range(kd,la,dist_cut);

	/* if there are some pairs, then tell us about them */
	if (kd_res_size(res)>0) {
	  while( !kd_res_end( res ) && matching_pairs<max_matches) {
	    /* get the data and position of the current result item */
	    data = (int*)kd_res_item( res, pos );
	    if (verbose>0) { 
	      diff=hypot(pos[0]-la[0],pos[1]-la[1]);
	      printf("# diff= %g\n",diff); 
	    }
	    if (verbose>2) {
	      printf("xp1[%d]= %g, xp2[%d]= %g xp1[%d]= %g, xp2[%d]= %g %g\n",
		     data[0],xp1[data[0]],datahold[0],xp2[datahold[0]],
		     data[1],xp1[data[1]],datahold[1],xp2[datahold[1]],la[0]);
	    }
	    pairoutput(data[0],data[1],datahold[0],datahold[1]);
	    /* go to the next entry */
	    kd_res_next( res );
	  }
	}
	/* free the results structure */
	kd_res_free(res);
	if (matching_pairs>max_matches) { i=j=n2; }
    }
  }

  /* output the best cooefficients */
  if (nbest>0) {
    printf("Transformation that fits the most asterisms:\n-t %g 0 %g 0 %g %g\n",
	   x1_factor,bestcoeff[0],y1_factor,bestcoeff[1]);
  } else {
    printf("No transformations with multiple asterisms found.\n");
  }

  /* free the tree */
  kd_free(kd);
  kd_free(kd_good);
  free ( (void *) fs1);
  free ( (void *) fs2);
  free ((void *) xp1);
  free ((void *) yp1);
  free ((void *) xp2);
  free ((void *) yp2);
  return 0;
}
