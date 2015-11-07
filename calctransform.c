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
#include "calctransform.h"

void
calctransform(double xp1[], double yp1[], double res[], int index1[], int index2[], int npoint, double param[])
{
  double s01, s02, s03, s11, s12, s13, s22, s23, s33;
  double d;
  int i;

  /* if there are fewer than three points, than the determinant will be zero. */
  if (npoint<3) {
    return;
  }

  s01=s02=s03=s11=s12=s13=s22=s23=s33=0;
  for (i=0;i<npoint;i++) {
    s01+=xp1[index1[i]]*res[index2[i]];
    s02+=yp1[index1[i]]*res[index2[i]];
    s03+=res[index2[i]];
    s11+=xp1[index1[i]]*xp1[index1[i]];
    s12+=xp1[index1[i]]*yp1[index1[i]];
    s13+=xp1[index1[i]];
    s22+=yp1[index1[i]]*yp1[index1[i]];
    s23+=yp1[index1[i]];
    s33++;
  }
  d=(s13*s13*s22-2*s12*s13*s23+s11*s23*s23+s12*s12*s33-s11*s22*s33);
  param[0]=((s03*s13*s22 - s03*s12*s23 - s02*s13*s23 + s01*s23*s23  + s02*s12*s33 - s01*s22*s33)/d);
  param[1]=((-s03*s12*s13 + s02*s13*s13  + s03*s11*s23 - s01*s13*s23 - s02*s11*s33 + s01*s12*s33)/d);
  param[2]=((s03*s12*s12  - s02*s12*s13 - s03*s11*s22 + s01*s13*s22 + s02*s11*s23 - s01*s12*s23)/d);
}

void
calctransform_noindex(double xp1[], double yp1[], double res[], int npoint, double param[])
{
  double s01, s02, s03, s11, s12, s13, s22, s23, s33;
  double d;
  int i;


  /* if there are fewer than three points, than the determinant will be zero. */
  if (npoint<3) {
    return;
  }

  s01=s02=s03=s11=s12=s13=s22=s23=s33=0;
  for (i=0;i<npoint;i++) {
    s01+=xp1[i]*res[i];
    s02+=yp1[i]*res[i];
    s03+=res[i];
    s11+=xp1[i]*xp1[i];
    s12+=xp1[i]*yp1[i];
    s13+=xp1[i];
    s22+=yp1[i]*yp1[i];
    s23+=yp1[i];
    s33++;
  }
  d=(s13*s13*s22-2*s12*s13*s23+s11*s23*s23+s12*s12*s33-s11*s22*s33);
  param[0]=((s03*s13*s22 - s03*s12*s23 - s02*s13*s23 + s01*s23*s23  + s02*s12*s33 - s01*s22*s33)/d);
  param[1]=((-s03*s12*s13 + s02*s13*s13  + s03*s11*s23 - s01*s13*s23 - s02*s11*s33 + s01*s12*s33)/d);
  param[2]=((s03*s12*s12  - s02*s12*s13 - s03*s11*s22 + s01*s13*s22 + s02*s11*s23 - s01*s12*s23)/d);
}
