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
#include "kdtree.h"
#include "loadfile.h"
#include "calctransform.h"


unsigned int n1, n2;
double *xp1, *yp1, *xp2, *yp2;
int listswapped=0, verbose=0, noswap=0, matching_pairs;
struct kdtree *kd_good;
double dist_cut=1e-5, trans_cut=1e-3, param2_factor=1000;
double bestcoeff[6];
int nbest=0;

typedef 
struct pair_data {
  double length;
  int i, j, k;
} pdata;

int
dcomp(const void *a, const void *b) {
  if (*((double *) a) <*((double *) b)) return 1;
  if (*((double *) a) >*((double *) b)) return -1;
  return 0;
}

int
pcomp(const void *a, const void *b) {
  if ( ((pdata *) a)->length <((pdata *)b)->length) return 1;
if ( ((pdata *) a)->length >((pdata *)b)->length) return -1;
  return 0;
}


int
addtriangle(double *param, int *tridata) {
  struct kdres *res;
  static int triangle_added;
  double pos[3];
  int *data, i;

  matching_pairs=0;
  if (triangle_added) {
    res=kd_nearest_range3(kd_good,param[0],param[1],param[2]/param2_factor,trans_cut);
    /* if there are some triangles with matching transforms, then tell us about them */
    if (kd_res_size(res)>0) {
      double sx1, sx2, sx3, sy1, sy2, sy3, s11, s12, s13, s22, s23, s33, d, coeff[6];
      sy1=sy2=sy3=sx1=sx2=sx3=s11=s12=s13=s22=s23=s33=0;

      for (i=0;i<3;i++) {
	if (listswapped) {
	  sx1+=xp2[tridata[i]]*xp1[tridata[i+3]];
	  sx2+=yp2[tridata[i]]*xp1[tridata[i+3]];
	  sx3+=xp1[tridata[i+3]];
	  sy1+=xp2[tridata[i]]*yp1[tridata[i+3]];
	  sy2+=yp2[tridata[i]]*yp1[tridata[i+3]];
	  sy3+=yp1[tridata[i+3]];
	  s11+=xp2[tridata[i]]*xp2[tridata[i]];
	  s12+=xp2[tridata[i]]*yp2[tridata[i]];
	  s13+=xp2[tridata[i]];
	  s22+=yp2[tridata[i]]*yp2[tridata[i]];
	  s23+=yp2[tridata[i]];
	} else {
	  sx1+=xp1[tridata[i]]*xp2[tridata[i+3]];
	  sx2+=yp1[tridata[i]]*xp2[tridata[i+3]];
	  sx3+=xp2[tridata[i+3]];
	  sy1+=xp1[tridata[i]]*yp2[tridata[i+3]];
	  sy2+=yp1[tridata[i]]*yp2[tridata[i+3]];
	  sy3+=yp2[tridata[i+3]];
	  s11+=xp1[tridata[i]]*xp1[tridata[i]];
	  s12+=xp1[tridata[i]]*yp1[tridata[i]];
	  s13+=xp1[tridata[i]];
	  s22+=yp1[tridata[i]]*yp1[tridata[i]];
	  s23+=yp1[tridata[i]];
	}
	s33++;
      }

      matching_pairs=kd_res_size(res);
      if (verbose>=0) {
	printf("x-transform: x2= %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
	printf("Number of matching triangle pairs: %d\n",matching_pairs);
      } 
      while( !kd_res_end( res ) ) {
	/* get the data and position of the current result item */
	data = (int*)kd_res_item( res, pos );
	if (verbose>=0) {
	  printf("Triangle pair with matching x-transform: x2= %g x1 + %g y1 + %g: {",pos[0],pos[1],pos[2]*param2_factor);
	  for (i=0;i<3;i++) {
	    printf(" %d",data[i]);
	  }
	  printf("} -> {");
	  for (i=3;i<6;i++) {
	    printf(" %d",data[i]);
	  }
	  printf("}\n");
	}
	for (i=0;i<3;i++) {
	  if (listswapped) {
	    sx1+=xp2[data[i]]*xp1[data[i+3]];
	    sx2+=yp2[data[i]]*xp1[data[i+3]];
	    sx3+=xp1[data[i+3]];
	    sy1+=xp2[data[i]]*yp1[data[i+3]];
	    sy2+=yp2[data[i]]*yp1[data[i+3]];
	    sy3+=yp1[data[i+3]];
	    s11+=xp2[data[i]]*xp2[data[i]];
	    s12+=xp2[data[i]]*yp2[data[i]];
	    s13+=xp2[data[i]];
	    s22+=yp2[data[i]]*yp2[data[i]];
	    s23+=yp2[data[i]];
	  } else {
	    sx1+=xp1[data[i]]*xp2[data[i+3]];
	    sx2+=yp1[data[i]]*xp2[data[i+3]];
	    sx3+=xp2[data[i+3]];
	    sy1+=xp1[data[i]]*yp2[data[i+3]];
	    sy2+=yp1[data[i]]*yp2[data[i+3]];
	    sy3+=yp2[data[i+3]];
	    s11+=xp1[data[i]]*xp1[data[i]];
	    s12+=xp1[data[i]]*yp1[data[i]];
	    s13+=xp1[data[i]];
	    s22+=yp1[data[i]]*yp1[data[i]];
	    s23+=yp1[data[i]];
	  }
	  s33++;
	}
	/* go to the next entry */
	kd_res_next( res );
      }
      d=(s13*s13*s22-2*s12*s13*s23+s11*s23*s23+s12*s12*s33-s11*s22*s33);
      coeff[0]=((sx3*s13*s22 - sx3*s12*s23 - sx2*s13*s23 + sx1*s23*s23  + sx2*s12*s33 - sx1*s22*s33)/d);
      coeff[1]=((-sx3*s12*s13 + sx2*s13*s13  + sx3*s11*s23 - sx1*s13*s23 - sx2*s11*s33 + sx1*s12*s33)/d);
      coeff[2]=((sx3*s12*s12  - sx2*s12*s13 - sx3*s11*s22 + sx1*s13*s22 + sx2*s11*s23 - sx1*s12*s23)/d);
      coeff[3]=((sy3*s13*s22 - sy3*s12*s23 - sy2*s13*s23 + sy1*s23*s23  + sy2*s12*s33 - sy1*s22*s33)/d);
      coeff[4]=((-sy3*s12*s13 + sy2*s13*s13  + sy3*s11*s23 - sy1*s13*s23 - sy2*s11*s33 + sy1*s12*s33)/d);
      coeff[5]=((sy3*s12*s12  - sy2*s12*s13 - sy3*s11*s22 + sy1*s13*s22 + sy2*s11*s23 - sy1*s12*s23)/d);
      if (verbose>=0) {
	printf("x2 = %g x1 + %g y1 + %g\n",
	       coeff[0],coeff[1],coeff[2]);
	printf("y2 = %g x1 + %g y1 + %g\n -t",
	       coeff[3],coeff[4],coeff[5]);
	for (i=0;i<6;i++) {
	  printf(" %g",coeff[i]);
	}
	printf("\n");
      }
      if (matching_pairs>nbest) {
	nbest=matching_pairs;
	for (i=0;i<6;i++) {
	  bestcoeff[i]=coeff[i];
	}
      }
    }
    /* free the results structure */
    kd_res_free(res);
  }
  triangle_added=1;

  return kd_insert3(kd_good, param[0], param[1], param[2]/param2_factor, (void *) tridata);
}

void
triangleoutput(int ih, int jh, int kh, int iah, int jah, int kah) {
  pdata pa[3], pb[3];
  int i;
  int index1[3], index2[3], *data;
  double param[3];

  /* calculate side lengths from short list and include points */
  pa[0].length=hypot(xp1[iah]-xp1[jah],yp1[iah]-yp1[jah]);
  pa[0].i=iah;
  pa[0].j=jah;
  pa[0].k=kah;
  pa[1].length=hypot(xp1[jah]-xp1[kah],yp1[jah]-yp1[kah]);
  pa[1].i=jah;
  pa[1].j=kah;
  pa[1].k=iah;
  pa[2].length=hypot(xp1[iah]-xp1[kah],yp1[iah]-yp1[kah]);
  pa[2].i=iah;
  pa[2].j=kah;
  pa[2].k=jah;
  /* sort them and include the point data */
  qsort((void *) pa,3,sizeof(pdata),pcomp);

  /* calculate side lengths from short list and include points */
  pb[0].length=hypot(xp2[ih]-xp2[jh],yp2[ih]-yp2[jh]);
  pb[0].i=ih;
  pb[0].j=jh;
  pb[0].k=kh;
  pb[1].length=hypot(xp2[jh]-xp2[kh],yp2[jh]-yp2[kh]);
  pb[1].i=jh;
  pb[1].j=kh;
  pb[1].k=ih;
  pb[2].length=hypot(xp2[ih]-xp2[kh],yp2[ih]-yp2[kh]);
  pb[2].i=ih;
  pb[2].j=kh;
  pb[2].k=jh;
  /* sort them and include the point data */
  qsort((void *) pb,3,sizeof(pdata),pcomp);

  /* determine point correspondences */
  for (i=0;i<3;i++) {
    index1[i]=pa[i].k;
    index2[i]=pb[i].k;
  }

  /* output data including the possibility that the lists were swapped */

  if (listswapped) { 
    /* output the pair information */
    if (verbose>1) {
      for (i=0;i<3;i++) {
	printf("Pair #%d %3d - %3d %10.4f ; %3d - %3d %10.4f\n",i,pb[i].i,pb[i].j,pb[i].length,pa[i].i,pa[i].j,pa[i].length);
      }
    }

    /* output the point correspondances */
    if (verbose>0) {
      printf("{%d %d %d} -> {%d %d %d}\n",index2[0],index2[1],index2[2],index1[0],index1[1],index1[2]);
    }

    /* calculate forward transformation */
    calctransform(xp2,yp2,xp1,index2,index1,3,param);
    data=(int *) malloc(sizeof(int)*6);
    for (i=0;i<3;i++) {
      data[i]=index2[i]; 
      data[i+3]=index1[i];
    }

    /* add it to the tree */
    if (addtriangle(param,data)) {
      printf("Unable to insert triangle into the tree at  %s:%d\n",__FILE__,__LINE__);
      return;
    }

    if (verbose) {
      printf("Transforms from 1->2\n");
      printf("x2 = %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
      calctransform(xp2,yp2,yp1,index2,index1,3,param);
      printf("y2 = %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
    }
    
    /* output reverse transformation */
    if (verbose>1) {
      printf("Transforms from 2->1\n");
      calctransform(xp1,yp1,xp2,index1,index2,3,param);
      printf("x1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
      calctransform(xp1,yp1,yp2,index1,index2,3,param);
      printf("y1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
    }
  } else {
    /* output the pair information */
    if (verbose>1) {
      for (i=0;i<3;i++) {
	printf("Pair #%d %3d - %3d %10.4f ; %3d - %3d %10.4f\n",i,pa[i].i,pa[i].j,pa[i].length,pb[i].i,pb[i].j,pb[i].length);
      }
    }
    /* output the point correspondances */
    if (verbose>0) {
      printf("{%d %d %d} -> {%d %d %d}\n",index1[0],index1[1],index1[2],index2[0],index2[1],index2[2]);
    }

    /* calculate forward transformation */
    calctransform(xp1,yp1,xp2,index1,index2,3,param);
    data=(int *) malloc(sizeof(int)*6);
    for (i=0;i<3;i++) {
      data[i]=index1[i]; 
      data[i+3]=index2[i];
    }

    /* add it to the tree */
    if (addtriangle(param,data)) {
      printf("Unable to insert triangle into the tree at  %s:%d\n",__FILE__,__LINE__);
      return;
    }

    if (verbose>0) {
      printf("Transforms from 1->2\n");
      printf("x2 = %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
      calctransform(xp1,yp1,yp2,index1,index2,3,param);
      printf("y2 = %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
    }
    
    /* output reverse transformation */
    if (verbose>1) {
      printf("Transforms from 2->1\n");
      calctransform(xp2,yp2,xp1,index2,index1,3,param);
      printf("x1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
      calctransform(xp2,yp2,yp1,index2,index1,3,param);
      printf("y1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
    }
  }  
}



int
main(int argc, char *argv[])
{
  int i, j, k, max_matches=20;
  int *data;
  double la[3], diff, ratioarray[2], pos[2],  atof();
  struct kdtree *kd;
  struct kdres *res;
  unsigned int cols1[]={1,2}, cols2[]={1,2};
  char **argptr, *filename1=NULL, *filename2=NULL;
  double *dptr[2];
  char *fs1, *fs2;

  fs1 = strdup(" \t");
  fs2 = strdup(" \t");

  if (argc<3) {
    printf("Format:\n\n   triangle_kd file1 file2 [options]\n\n\
   where the options can appear anywhere in any order:\n\n\
   -d distance_cutoff  how small of a distance to call a match - default %g\n\
   -t transform_cutoff how small of a distance to call a match - default %g\n\
   -p translate_factor factor to scale the x-translation       - default %g\n\
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
   triangle_kd will try to match triangles in the two catalogues.  It is useful\n\
   when one expects there to be rotation, translation and scaling between the\n\
   two catalogues but no shearing.  If you do expect shearing, use quad_kd.\n",
	   dist_cut,trans_cut,param2_factor,max_matches,cols1[0],cols1[1],cols2[0],cols2[1]);
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
    } else if (strstr(*argptr,"-p")) {
      if (++argptr<argv+argc) {
	param2_factor=atof(*argptr);
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


  /* output some information at the top */
  if (verbose>0) {
    printf("# List 1: %s has %u points and %lu triangles\n# List 2: %s has %u points and %lu triangles\n",
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


  /* create the kd-tree for the triangles*/
  kd = kd_create(2);
  /* designate a function to deallocate the data */
  kd_data_destructor(kd,free);

  /* build the tree */
  for (i=0;i<n1-2;i++) {
    for (j=i+1;j<n1-1;j++) {
      for (k=j+1;k<n1;k++) {
	/* calculate side lengths */
	la[0]=hypot(xp1[i]-xp1[j],yp1[i]-yp1[j]);
	la[1]=hypot(xp1[j]-xp1[k],yp1[j]-yp1[k]);
	la[2]=hypot(xp1[i]-xp1[k],yp1[i]-yp1[k]);
	/* sort them */
	qsort((void *) la,3,sizeof(double),dcomp);
	/* if the smallest side is less than a pixel, skip this triangle */
	if (la[2]<1) break;
	/* calculate the side ratio */
	ratioarray[0]=la[1]/la[0];
	ratioarray[1]=la[2]/la[0];
	/* allocate an array to hold the points */
	data=(int *) malloc(sizeof(int)*3);
	data[0]=i; data[1]=j; data[2]=k;
	/* add it to the tree */
	if (kd_insert(kd, ratioarray, (void *) data)) {
	  printf("Unable to insert point into the tree at %s:%d\n",__FILE__,__LINE__);
	  return -1;
	}
      }
    }
  }


  /* go through the triangles from the longer list */

  /* create the kd-tree for the matches */
  kd_good = kd_create(3);
  /* designate a function to deallocate the data */
  kd_data_destructor(kd_good,free);

  for (i=0;i<n2-2;i++) {
    for (j=i+1;j<n2-1;j++) {
      for (k=j+1;k<n2;k++) {
	/* calculate side lengths */
	la[0]=hypot(xp2[i]-xp2[j],yp2[i]-yp2[j]);
	la[1]=hypot(xp2[j]-xp2[k],yp2[j]-yp2[k]);
	la[2]=hypot(xp2[i]-xp2[k],yp2[i]-yp2[k]);
	/* sort them */
	qsort((void *) la,3,sizeof(double),dcomp);

	/* if the smallest side is less than a pixel, skip this triangle */
	if (la[2]<1) break;

	/* calculate the ratio */
	ratioarray[0]=la[1]/la[0];
	ratioarray[1]=la[2]/la[0];

	/* find all the triangles from the first (short) list that are within the dist_cut */
	res=kd_nearest_range(kd,ratioarray,dist_cut);

	/* if there are some triangles, then tell us about them */
	if (kd_res_size(res)>0) {
	  while( !kd_res_end( res ) && matching_pairs<max_matches) {
	    /* get the data and position of the current result item */
	    data = (int*)kd_res_item( res, pos );
	    if (verbose>0) { 
	      diff=hypot(pos[0]-ratioarray[0],pos[1]-ratioarray[1]);
	      printf("# diff= %g\n",diff); 
	    }
	    triangleoutput(i,j,k,data[0],data[1],data[2]);
	    /* go to the next entry */
	    kd_res_next( res );
	  }
	}
	/* free the results structure */
	kd_res_free(res);
	if (matching_pairs>max_matches) { i=j=k=n2; }
      }
    }
  }

  /* output the best cooefficients */
  if (nbest>0) {
    printf("Transformation that fits the most asterisms:\n-t");
    for (i=0;i<6;i++) {
      printf(" %g",bestcoeff[i]);
    }
    printf("\n");
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
