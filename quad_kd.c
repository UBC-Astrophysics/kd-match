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


int n1, n2;
double *xp1, *yp1, *xp2, *yp2;
int listswapped=0, verbose=0, noswap=0, matching_pairs;
double dist_cut=3e-3, trans_cut=1e-3, param2_factor=1000;
struct kdtree *kd_good;
double bestcoeff[6];
int nbest=0;

typedef 
struct quad_data {
  double area;
  int i, j, k, l;
} qdata;

int
dcomp(const void *a, const void *b) {
  if (*((double *) a) <*((double *) b)) return 1;
  if (*((double *) a) >*((double *) b)) return -1;
  return 0;
}

int
tcomp(const void *a, const void *b) {
  if ( ((qdata *) a)->area <((qdata *)b)->area) return 1;
if ( ((qdata *) a)->area >((qdata *)b)->area) return -1;
  return 0;
}


int
addquad(double *param, int *quaddata) {
  struct kdres *res;
  static int quad_added;
  double pos[3];
  int *data, i, quad_cnt;
  
  matching_pairs=0;
  if (quad_added) {
    res=kd_nearest_range3(kd_good,param[0],param[1],param[2]/param2_factor,trans_cut);
    /* if there are some quads with matching transforms, then tell us about them */
    if (kd_res_size(res)>0) {
      double sx1, sx2, sx3, sy1, sy2, sy3, s11, s12, s13, s22, s23, s33, d, coeff[6];
      sy1=sy2=sy3=sx1=sx2=sx3=s11=s12=s13=s22=s23=s33=0;

      for (i=0;i<4;i++) {
	if (listswapped) {
	  sx1+=xp2[quaddata[i]]*xp1[quaddata[i+4]];
	  sx2+=yp2[quaddata[i]]*xp1[quaddata[i+4]];
	  sx3+=xp1[quaddata[i+4]];
	  sy1+=xp2[quaddata[i]]*yp1[quaddata[i+4]];
	  sy2+=yp2[quaddata[i]]*yp1[quaddata[i+4]];
	  sy3+=yp1[quaddata[i+4]];
	  s11+=xp2[quaddata[i]]*xp2[quaddata[i]];
	  s12+=xp2[quaddata[i]]*yp2[quaddata[i]];
	  s13+=xp2[quaddata[i]];
	  s22+=yp2[quaddata[i]]*yp2[quaddata[i]];
	  s23+=yp2[quaddata[i]];
	} else {
	  sx1+=xp1[quaddata[i]]*xp2[quaddata[i+4]];
	  sx2+=yp1[quaddata[i]]*xp2[quaddata[i+4]];
	  sx3+=xp2[quaddata[i+4]];
	  sy1+=xp1[quaddata[i]]*yp2[quaddata[i+4]];
	  sy2+=yp1[quaddata[i]]*yp2[quaddata[i+4]];
	  sy3+=yp2[quaddata[i+4]];
	  s11+=xp1[quaddata[i]]*xp1[quaddata[i]];
	  s12+=xp1[quaddata[i]]*yp1[quaddata[i]];
	  s13+=xp1[quaddata[i]];
	  s22+=yp1[quaddata[i]]*yp1[quaddata[i]];
	  s23+=yp1[quaddata[i]];
	}
	s33++;
      }
      
      matching_pairs=kd_res_size(res);
      if (verbose>=0) {
	printf("x-transform: x2= %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
	printf("Number of matching quad pairs: %d\n",matching_pairs);
      }
      while( !kd_res_end( res ) ) {
	/* get the data and position of the current result item */
	data = (int*)kd_res_item( res, pos );
	if (verbose>=0) {
	  printf("Quad pair with matching x-transform: x2= %g x1 + %g y1 + %g: {",pos[0],pos[1],pos[2]*param2_factor);
	  for (i=0;i<4;i++) {
	    printf(" %d",data[i]);
	  }
	  printf("} -> {");
	  for (i=4;i<8;i++) {
	    printf(" %d",data[i]);
	  }
	  printf("}\n");
	}
	for (i=0;i<4;i++) {
	  if (listswapped) {
	    sx1+=xp2[data[i]]*xp1[data[i+4]];
	    sx2+=yp2[data[i]]*xp1[data[i+4]];
	    sx3+=xp1[data[i+4]];
	    sy1+=xp2[data[i]]*yp1[data[i+4]];
	    sy2+=yp2[data[i]]*yp1[data[i+4]];
	    sy3+=yp1[data[i+4]];
	    s11+=xp2[data[i]]*xp2[data[i]];
	    s12+=xp2[data[i]]*yp2[data[i]];
	    s13+=xp2[data[i]];
	    s22+=yp2[data[i]]*yp2[data[i]];
	    s23+=yp2[data[i]];
	  } else {
	    sx1+=xp1[data[i]]*xp2[data[i+4]];
	    sx2+=yp1[data[i]]*xp2[data[i+4]];
	    sx3+=xp2[data[i+4]];
	    sy1+=xp1[data[i]]*yp2[data[i+4]];
	    sy2+=yp1[data[i]]*yp2[data[i+4]];
	    sy3+=yp2[data[i+4]];
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
	printf("x2 = %g x1 + %g y1 + %g\n",coeff[0],coeff[1],coeff[2]);
	printf("y2 = %g x1 + %g y1 + %g\n -t",coeff[3],coeff[4],coeff[5]);
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
  quad_added=1;

  return kd_insert3(kd_good, param[0], param[1], param[2]/param2_factor, (void *) quaddata);
}

void
quadoutput(int ih, int jh, int kh, int lh, 
	   int iah, int jah, int kah, int lah) {
  qdata pa[4], pb[4];
  int i, index1[4], index2[4], *data;
  double param[3];

  /* calculate triangle areas from short list and include points */

  pa[0].area=fabs((xp1[iah]-xp1[jah])*(yp1[jah]-yp1[kah])-(xp1[jah]-xp1[kah])*(yp1[iah]-yp1[jah]));  /* leave out l */
  pa[0].i=iah;
  pa[0].j=jah;
  pa[0].k=kah;
  pa[0].l=lah;

  pa[1].area=fabs((xp1[iah]-xp1[jah])*(yp1[jah]-yp1[lah])-(xp1[jah]-xp1[lah])*(yp1[iah]-yp1[jah]));  /* leave out k */
  pa[1].i=iah;
  pa[1].j=jah;
  pa[1].k=lah;
  pa[1].l=kah;

  pa[2].area=fabs((xp1[iah]-xp1[lah])*(yp1[lah]-yp1[kah])-(xp1[lah]-xp1[kah])*(yp1[iah]-yp1[lah]));  /* leave out j */
  pa[2].i=iah;
  pa[2].j=kah;
  pa[2].k=lah;
  pa[2].l=jah;

  pa[3].area=fabs((xp1[lah]-xp1[jah])*(yp1[jah]-yp1[kah])-(xp1[jah]-xp1[kah])*(yp1[lah]-yp1[jah]));  /* leave out i */
  pa[3].i=jah;
  pa[3].j=kah;
  pa[3].k=lah;
  pa[3].l=iah;

  /* sort them and include the point data */
  qsort((void *) pa,4,sizeof(qdata),tcomp);

  /* calculate triangle areas from long list and include points */

  pb[0].area=fabs((xp2[ih]-xp2[jh])*(yp2[jh]-yp2[kh])-(xp2[jh]-xp2[kh])*(yp2[ih]-yp2[jh]));  /* leave out l */
  pb[0].i=ih;
  pb[0].j=jh;
  pb[0].k=kh;
  pb[0].l=lh;

  pb[1].area=fabs((xp2[ih]-xp2[jh])*(yp2[jh]-yp2[lh])-(xp2[jh]-xp2[lh])*(yp2[ih]-yp2[jh]));  /* leave out k */
  pb[1].i=ih;
  pb[1].j=jh;
  pb[1].k=lh;
  pb[1].l=kh;

  pb[2].area=fabs((xp2[ih]-xp2[lh])*(yp2[lh]-yp2[kh])-(xp2[lh]-xp2[kh])*(yp2[ih]-yp2[lh]));  /* leave out j */
  pb[2].i=ih;
  pb[2].j=kh;
  pb[2].k=lh;
  pb[2].l=jh;

  pb[3].area=fabs((xp2[lh]-xp2[jh])*(yp2[jh]-yp2[kh])-(xp2[jh]-xp2[kh])*(yp2[lh]-yp2[jh]));  /* leave out i */
  pb[3].i=jh;
  pb[3].j=kh;
  pb[3].k=lh;
  pb[3].l=ih;

  /* sort them and include the point data */
  qsort((void *) pb,4,sizeof(qdata),tcomp);

  /* determine point correspondences */

  for (i=0;i<4;i++) {
    index1[i]=pa[i].l;
    index2[i]=pb[i].l;
  }

  /* output data including the possibility that the lists were swapped */

  if (listswapped) { 
    /* output the triangle information */
    if (verbose>1) {
      for (i=0;i<4;i++) {
	printf("Tri #%d %3d - %3d - %3d %10.4f ; %3d - %3d - %3d %10.4f\n",i,pb[i].i,pb[i].j,pb[i].k,pb[i].area,pa[i].i,pa[i].j,pa[i].k,pa[i].area);
      }
    }
    /* output the point correspondances */
    if (verbose>0) {
      printf("{%d %d %d %d} -> {%d %d %d %d}\n",index2[0],index2[1],index2[2],index2[3],index1[0],index1[1],index1[2],index1[3]);
    }

    /* calculate forward transformation */
    calctransform(xp2,yp2,xp1,index2,index1,4,param);
    data=(int *) malloc(sizeof(int)*8);
    for (i=0;i<4;i++) {
      data[i]=index2[i]; 
      data[i+4]=index1[i];
    }

    /* add it to the tree */
    if (addquad(param,data)) {
      printf("Unable to insert quad into the tree at  __FILE__:__LINENO__\n");
      return;
    }

    if (verbose>0) {
      printf("Transforms from 1->2\n");
      printf("x2 = %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
      calctransform(xp2,yp2,yp1,index2,index1,4,param);
      printf("y2 = %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
    }
    
    /* output reverse transformation */
    if (verbose>1) {
      printf("Transforms from 2->1\n");
      calctransform(xp1,yp1,xp2,index1,index2,4,param);
      printf("x1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
      calctransform(xp1,yp1,yp2,index1,index2,4,param);
      printf("y1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
    }
  } else {
    /* output the triangle information */
    if (verbose>1) {
      for (i=0;i<4;i++) {
	printf("Tri #%d %3d - %3d - %3d %10.4f ; %3d - %3d - %3d %10.4f\n",i,pa[i].i,pa[i].j,pa[i].k,pa[i].area,pb[i].i,pb[i].j,pb[i].k,pb[i].area);
      }
    }

    /* output the point correspondances */
    if (verbose>0) {
      printf("{%d %d %d %d} -> {%d %d %d %d}\n",index1[0],index1[1],index1[2],index1[3],index2[0],index2[1],index2[2],index2[3]);
    }

    /* calculate forward transformation */
    calctransform(xp1,yp1,xp2,index1,index2,4,param);
    data=(int *) malloc(sizeof(int)*8);
    for (i=0;i<4;i++) {
      data[i]=index1[i]; 
      data[i+4]=index2[i];
    }

    /* add it to the tree */
    if (addquad(param,data)) {
      printf("Unable to insert quad into the tree at  __FILE__:__LINENO__\n");
      return;
    }

    if (verbose>0) {
      printf("Transforms from 1->2\n");
      printf("x2 = %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
      calctransform(xp1,yp1,yp2,index1,index2,4,param);
      printf("y2 = %g x1 + %g y1 + %g\n",param[0],param[1],param[2]);
    }
    
    /* output reverse transformation */
    if (verbose>1) {
      printf("Transforms from 2->1\n");
      calctransform(xp2,yp2,xp1,index2,index1,4,param);
      printf("x1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
      calctransform(xp2,yp2,yp1,index2,index1,4,param);
      printf("y1 = %g x2 + %g y2 + %g\n",param[0],param[1],param[2]);
    }
  }
}

int
main(int argc, char *argv[])
{
  FILE *in;
  int i, j, k, l, ih, jh, kh, lh;
  int iah, jah, kah, lah, *data, max_matches=20;
  double la[4], bestdiff, diff, ratioarray[2], pos[2], atof();
  unsigned int cols1[]={1,2}, cols2[]={1,2};
  char **argptr, *filename1=NULL, *filename2=NULL;
  struct kdtree *kd;
  struct kdres *res;
  double *dptr[2];
  char *fs1, *fs2;

  fs1 = strdup(" \t");
  fs2 = strdup(" \t");

  if (argc<3) {
    printf("Format:\n\n   quad_kd file1 file2 [options]\n\n\
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
   quad_kd will try to match quadrilaterals in the two catalogues.  It is useful\n\
   when one expects there to be shearing as well as rotation, translation and\n\
   scaling between the two catalogues.  If you don't expect shearing, use\n\
   triangle_kd.\n",
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
    printf("# List 1: %s has %d points and %ld quads\n# List 2: %s has %d points and %ld quads\n# dist_cut= %g\n%s",
	 argv[1],n1,n1*(n1-1l)*(n1-2l)*(n1-3l)/24,argv[2],n2,n2*(n2-1l)*(n2-2l)*(n1-3l)/24,dist_cut,(listswapped ? "# For internal use the lists were swapped.\n ": ""));
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
  } else {
    listswapped=0;
  }

  /* create the kd-tree */
  kd = kd_create(2);
  /* designate a function to deallocate the data */
  kd_data_destructor(kd,free);

  /* build the tree */
  for (i=0;i<n1-3;i++) {
    for (j=i+1;j<n1-2;j++) {
      for (k=j+1;k<n1-1;k++) {
	for (l=k+1;l<n1;l++) {
	/* calculate triangle areas */
	la[0]=fabs((xp1[i]-xp1[j])*(yp1[j]-yp1[k])-(xp1[j]-xp1[k])*(yp1[i]-yp1[j]));  /* leave out l */
	la[1]=fabs((xp1[i]-xp1[j])*(yp1[j]-yp1[l])-(xp1[j]-xp1[l])*(yp1[i]-yp1[j]));  /* leave out k */
	la[2]=fabs((xp1[i]-xp1[l])*(yp1[l]-yp1[k])-(xp1[l]-xp1[k])*(yp1[i]-yp1[l]));  /* leave out j */
	la[3]=fabs((xp1[l]-xp1[j])*(yp1[j]-yp1[k])-(xp1[j]-xp1[k])*(yp1[l]-yp1[j]));  /* leave out i */
	/* sort them */
	qsort((void *) la,4,sizeof(double),dcomp);
	/* if the smallest area is less than a pixel, skip this quad */
	if (la[3]<1) break;
	/* calculate the area ratio for the 2nd biggest to biggest */
        /* and 3rd biggest to biggest, the final ratio is determined */
	/* by the other two */
	ratioarray[0]=la[1]/la[0];
	ratioarray[1]=la[2]/la[0];
	/* allocate an array to hold the points */
	if ( (data=(int *) malloc(sizeof(int)*4))==NULL) {
	  printf("Unable to allocate data at __FILE__:__LINENO__\n");
	  return -1;
	}
	data[0]=i; data[1]=j; data[2]=k; data[3]=l;
	/* add it to the tree */
	if (kd_insert(kd, ratioarray, (void *) data)) {
	  printf("Unable to insert ratio into the tree at  __FILE__:__LINENO__\n");
	  return -1;

	}
      }
      }
    }
  }


  /* go through the quads from the longer list */

  /* create the kd-tree for the matches */
  kd_good = kd_create(3);
  /* designate a function to deallocate the data */
  kd_data_destructor(kd_good,free);

  for (i=0;i<n2-3;i++) {
    for (j=i+1;j<n2-2;j++) {
      for (k=j+1;k<n2-1;k++) {
	for (l=k+1;l<n2;l++) {
	  /* calculate triangle areas */
	  la[0]=fabs((xp2[i]-xp2[j])*(yp2[j]-yp2[k])-(xp2[j]-xp2[k])*(yp2[i]-yp2[j]));  /* leave out l */
	  la[1]=fabs((xp2[i]-xp2[j])*(yp2[j]-yp2[l])-(xp2[j]-xp2[l])*(yp2[i]-yp2[j]));  /* leave out k */
	  la[2]=fabs((xp2[i]-xp2[l])*(yp2[l]-yp2[k])-(xp2[l]-xp2[k])*(yp2[i]-yp2[l]));  /* leave out j */
	  la[3]=fabs((xp2[l]-xp2[j])*(yp2[j]-yp2[k])-(xp2[j]-xp2[k])*(yp2[l]-yp2[j]));  /* leave out i */
	  /* sort them */
	  qsort((void *) la,4,sizeof(double),dcomp);
	  
	  /* if the smallest area is less than a pixel, skip this quad */
	  if (la[3]<1) break;
	  /* calculate the area ratio for the 2nd biggest to biggest */
	  /* and 3rd biggest to biggest, the final ratio is determined */
	  /* by the other two */
	  ratioarray[0]=la[1]/la[0];
	  ratioarray[1]=la[2]/la[0];
	  
	  /* find all the quads from the first (short) list that are within the dist_cut */
	  res=kd_nearest_range(kd,ratioarray,dist_cut);
	  
	  /* if there are some quads, then tell us about them */
	  if (kd_res_size(res)>0) {
	    while( !kd_res_end( res ) && matching_pairs < max_matches) {
	      /* get the data and position of the current result item */
	      data = (int*) kd_res_item( res, pos );
	      if (verbose>0) {
		diff=hypot(pos[0]-ratioarray[0],pos[1]-ratioarray[1]);
		printf("# %g %g\n",ratioarray[0],ratioarray[1]);
		printf("# %g %g\n",pos[0],pos[1]);
		printf("# diff= %g  %d %d %d %d %d %d %d %d\n",diff,i,j,k,l,data[0],data[1],data[2],data[3]);
	      }
	      quadoutput(i,j,k,l,data[0],data[1],data[2],data[3]);
	      /* go to the next entry */
	      kd_res_next( res );
	    }
	  }
	  /* free the results structure */
	  kd_res_free(res);
	  if (matching_pairs>max_matches) { i=j=k=l=n2; }
	}
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
