#include <tconcoord.h>


/*===================================================*/
real dist2(rvec v1, rvec v2)
/* return squared distance */
{
  rvec diff;
  rvec_sub(v1,v2,diff);
  return norm2(diff);
}
/*===================================================*/
real dist_ij(rvec v1, rvec v2)
/* return distance */
{
  rvec diff;
  rvec_sub(v1,v2,diff);
  return norm(diff);
}
/*===================================================*/

/* real torsion(rvec v1,rvec v2, rvec v3, rvec v4) */

/* { */
/*   rvec ij,jk,kl; */
/*   rvec ij_jk,jk_kl; */
/*   real phi; */
/*   rvec_sub(v2,v1,ij); */
/*   rvec_sub(v3,v2,jk); */
/*   rvec_sub(v4,v3,kl); */
/*   vectorProduct(ij,jk,ij_jk); */
/*   vectorProduct(jk,kl,jk_kl); */
/*   phi=angle_ij(ij_jk,jk_kl); */
/*   if(determinant(ij_jk,jk_kl,jk)<=0.) */
/*     phi=2.*PI-phi; */
/*   return phi; */
/* } */
/*===================================================*/
real dihedral(const rvec xi,const rvec xj,const rvec xk,const rvec xl)
/* gromacs dihedral angle */          
{
  real ipr,phi,cos_phi,sign;
  rvec r_ij,r_kj,r_kl,m,n;
  
  rvec_sub(xi,xj,r_ij);  
  rvec_sub(xj,xk,r_kj);  
  rvec_sub(xk,xl,r_kl);  

  oprod(r_ij,r_kj,m);        
  oprod(r_kj,r_kl,n);        
  cos_phi=cos_angle(m,n);
  if(cos_phi < -1.) cos_phi = -1.;
  else if(cos_phi > 1.) cos_phi = 1.;
  phi=acos(cos_phi);         
  ipr=iprod(r_ij,n);         
  sign=(ipr<0.0)?-1.0:1.0;
  phi=sign*phi;              
                             
  return phi;
}
/*===================================================*/
real angle_ij_ik(rvec v1,rvec v2, rvec v3)
{
  real cos_phi;
  rvec ij,ik;
  real absij,absik;
  vectorMinusVector(v1,v2,ij);
  vectorMinusVector(v1,v3,ik);
  absij=absVec(ij);
  absik=absVec(ik);
  cos_phi=scalarProduct(ij,ik)/(absij*absik);
  if(cos_phi < -1.) cos_phi = -1.;
  else if(cos_phi > 1.) cos_phi = 1.;
  return acos(cos_phi);
}
/*===============================================================*/
real scalarProduct(rvec v1, rvec v2)
{
  real prod;
  
  prod = (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
  return prod;
}
/*===============================================================*/
real absVec(rvec v)
{
  real x;
  x = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  return x;
}
/*===============================================================*/
void vectorProduct(rvec v1, rvec v2, rvec v)
{
  v[0] = v1[1]*v2[2]-v1[2]*v2[1];
  v[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v[2] = v1[0]*v2[1]-v1[1]*v2[0];
}
/*===============================================================*/
void vectorPlusVector(rvec v1,rvec v2,rvec v)
{
  v[0] = v1[0]  + v2[0];
  v[1] = v1[1]  + v2[1];
  v[2] = v1[2]  + v2[2];
}
/*===============================================================*/
void vectorMinusVector(rvec v1,rvec v2,rvec v)
{
  v[0] = v1[0]  - v2[0];
  v[1] = v1[1]  - v2[1];
  v[2] = v1[2]  - v2[2];
}
/*===============================================================*/
real determinant(rvec a,rvec b, rvec c)
{
  real det;
  det = a[0]*(b[1]*c[2]-c[1]*b[2]) +
    b[0]*(c[1]*a[2]-a[1]*c[2]) +
    c[0]*(a[1]*b[2]-b[1]*a[2]);
  return det;
}

/*===============================================================*/

void scaleVector(rvec v_, rvec vsc_)
{
  real abs;
  
  abs = absVec(v_);
  
  vsc_[0] = v_[0]/abs;
  vsc_[1] = v_[1]/abs;
  vsc_[2] = v_[2]/abs;
  
}

/*===============================================================*/
void createRotationMatrix1(rvec a_,matrix m)
{
  m[0][0] = a_[0]*a_[0];
  m[0][1] = a_[0]*a_[1];
  m[0][2] = a_[0]*a_[2];
  
  m[1][0] = a_[1]*a_[0];
  m[1][1] = a_[1]*a_[1];
  m[1][2] = a_[1]*a_[2];
  
  m[2][0] = a_[2]*a_[0];
  m[2][1] = a_[2]*a_[1];
  m[2][2] = a_[2]*a_[2];
  
}
/*===============================================================*/
void createRotationMatrix2(rvec a_, matrix m)
{
  m[0][0] = 0.0;
  m[0][1] = -a_[2];
  m[0][2] = a_[1];
  
  m[1][0] = a_[2];
  m[1][1] = 0.0;
  m[1][2] = -a_[0];
  
  m[2][0] = -a_[1];
  m[2][1] = a_[0];
  m[2][2] = 0.0;

}
/*===============================================================*/

void MatrixVectorMultiply(matrix m,rvec x_, rvec v_)
{
  v_[0] = m[0][0]*x_[0] + m[0][1]*x_[1] + m[0][2]*x_[2];
  v_[1] = m[1][0]*x_[0] + m[1][1]*x_[1] + m[1][2]*x_[2];
  v_[2] = m[2][0]*x_[0] + m[2][1]*x_[1] + m[2][2]*x_[2];
}


/*===============================================================*/

void rotateVectorAroundVector(rvec v1_,rvec v2_,matrix tm1,
                              matrix tm2, real phi)
{
  rvec vec,a,b,c,d,e;

  vectorMinusVector(v1_,v2_,vec);
  MatrixVectorMultiply(tm1,vec,b);
  MatrixVectorMultiply(tm2,vec,d);
  
  a[0] = cos(phi) * vec[0];
  a[1] = cos(phi) * vec[1];
  a[2] = cos(phi) * vec[2];
  
  c[0] = -cos(phi) * b[0];
  c[1] = -cos(phi) * b[1];
  c[2] = -cos(phi) * b[2];
  
  e[0] = sin(phi) * d[0];
  e[1] = sin(phi) * d[1];
  e[2] = sin(phi) * d[2];
  
  vec[0] = a[0] + b[0] + c[0] + e[0];
  vec[1] = a[1] + b[1] + c[1] + e[1];
  vec[2] = a[2] + b[2] + c[2] + e[2];

  v1_[0] = v2_[0]+vec[0];
  v1_[1] = v2_[1]+vec[1];
  v1_[2] = v2_[2]+vec[2];
}  

/*===============================================================*/
void rotate_atom(rvec v1_, rvec v2_, rvec v3_, real phi)
{
  /* rotate vector 3 around v2-v1 */
  rvec rotvec,s_rotvec;
  vectorMinusVector(v2_,v1_,rotvec);
  scaleVector(rotvec,s_rotvec);
  matrix tm1,tm2;
  createRotationMatrix1(s_rotvec,tm1);
  createRotationMatrix2(s_rotvec,tm2);
  
  /* rotate vector 3 */
  rotateVectorAroundVector(v3_,v2_,tm1,tm2,phi);

}

/*===============================================================*/

void princ_comp(int n, rvec x[],
                    matrix trans,rvec d)
{

#define NDIM 4
  int  i,j,ai,m,nrot;
  real mm,rx,ry,rz;
  double **inten,dd[NDIM],tvec[NDIM],**ev;

  real temp;
  
  snew(inten,NDIM);
  snew(ev,NDIM);

  for(i=0; (i<NDIM); i++) {
    snew(inten[i],NDIM);
    snew(ev[i],NDIM);
    dd[i]=0.0;
  }
  
  for(i=0; (i<NDIM); i++)
    for(m=0; (m<NDIM); m++)
      inten[i][m]=0;
  for(i=0; (i<n); i++) {
    rx=x[i][XX];
    ry=x[i][YY];
    rz=x[i][ZZ];
    inten[0][0]+=(sqr(ry)+sqr(rz));
    inten[1][1]+=(sqr(rx)+sqr(rz));
    inten[2][2]+=(sqr(rx)+sqr(ry));
    inten[1][0]-=(ry*rx);
    inten[2][0]-=(rx*rz);
    inten[2][1]-=(rz*ry);
  }
  inten[0][1]=inten[1][0];
  inten[0][2]=inten[2][0];
  inten[1][2]=inten[2][1];
  
  for(i=0; (i<DIM); i++) {
    for(m=0; (m<DIM); m++)
      trans[i][m]=inten[i][m];
  }

  /* Call numerical recipe routines */
  jacobi(inten,3,dd,ev,&nrot);
  
  /* Sort eigenvalues in descending order */
#define SWAPPER(i) 			\
  if (fabs(dd[i+1]) > fabs(dd[i])) {	\
    temp=dd[i];			\
    for(j=0; (j<NDIM); j++) tvec[j]=ev[j][i];\
    dd[i]=dd[i+1];			\
    for(j=0; (j<NDIM); j++) ev[j][i]=ev[j][i+1];		\
    dd[i+1]=temp;			\
    for(j=0; (j<NDIM); j++) ev[j][i+1]=tvec[j];			\
  }
  SWAPPER(0)
  SWAPPER(1)
  SWAPPER(0)
    
  for(i=0; (i<DIM); i++) {
    d[i]=dd[i];
    for(m=0; (m<DIM); m++)
      trans[i][m]=ev[m][i]; 

  }
  for(i=0; (i<NDIM); i++) {
    sfree(inten[i]);
    sfree(ev[i]);
  }
  sfree(inten);
  sfree(ev);
}
/*===================================================*/
void rotate_atoms(int n,rvec x[],matrix trans)
{
  real   xt,yt,zt;
  int    i;
  
  for(i=0; (i<n); i++) {
    xt=x[i][XX];
    yt=x[i][YY];
    zt=x[i][ZZ];
    x[i][XX]=trans[XX][XX]*xt+trans[XX][YY]*yt+trans[XX][ZZ]*zt;
    x[i][YY]=trans[YY][XX]*xt+trans[YY][YY]*yt+trans[YY][ZZ]*zt;
    x[i][ZZ]=trans[ZZ][XX]*xt+trans[ZZ][YY]*yt+trans[ZZ][ZZ]*zt;
  }
}


/*===================================================*/

int findEulerAngles(matrix3x3 R, double *Rx, double *Ry, double *Rz)

{

	double theta1, theta2, psi1, psi2, phi1, phi2;


	if (FPNotEqualTo(R[2][0], 1) && FPNotEqualTo(R[2][0], -1)) {
		// Compute the two possible thetas
		theta1 = -asin(R[2][0]);
		theta2 = PI - theta1;

		// Compute the corresponding psi
		psi1 = atan2(R[2][1] / cos(theta1), R[2][2] / cos(theta1));
		psi2 = atan2(R[2][1] / cos(theta2), R[2][2] / cos(theta2));

		// Compute the corresponding phi
		phi1 = atan2(R[1][0] / cos(theta1), R[0][0] / cos(theta1));
		phi2 = atan2(R[1][0] / cos(theta2), R[0][0] / cos(theta2));


		// Set Euler angles to the first solution
		SetEulerAngles(psi1, theta1, phi1);
		return TRUE;
    } else {
		// Determine theta
      theta1 = 0;
      
		if (FPEqualTo(R[2][0], -1))
			theta1 = PI / 2;
		else if (FPEqualTo(R[2][0], 1))
			theta1 = -PI / 2;

		// Set phi = 0
		phi1 = 0;

		// Now compute psi
		psi1 = atan2(R[0][1], R[0][2]);
		SetEulerAngles(psi1, theta1, phi1);
		return TRUE;
		
	}
	return FALSE;
}

/*===================================================*/

real calc_rmsd(rvec *x1, rvec *x2, int size)
{
  int i;
  rvec diff;
  real rmsd = 0.;
  for(i=0;i<size;i++){
    rvec_sub(x1[i],x2[i],diff);
    rmsd+=norm(diff);
  }
  return (real) rmsd/size;
}

/*===================================================*/

void max_coords(t_atomlist *al, matrix box)
{
  clear_mat(box);
  int i,k;
  rvec max,min;
  for(i=0;i<DIM;i++){
    max[i] = -1000.;
    min[i] = 1000.;
  }
  for(i=0;i<al->natoms;i++){
    for(k=0;k<DIM;k++){
      if(al->x[i][k] < min[k]) {
        min[k] = al->x[i][k];
/*         *min = i; */
      }
      
      if(al->x[i][k] > max[k]) {
        max[k] = al->x[i][k];
/*         *max= i; */
      }
      
        
    }
  }
  for(i=0;i<DIM;i++){
    box[i][i] = max[i] - min[i];
  }
}

/*===================================================*/

void max_coords2(t_atomlist *al, matrix box, int *maxi, int *mini)
{
  clear_mat(box);
  int i,k;
  rvec max,min;
  for(i=0;i<DIM;i++){
    max[i] = -1000.;
    min[i] = 1000.;
  }
  for(i=0;i<al->natoms;i++){
    for(k=0;k<DIM;k++){
      if(al->x[i][k] < min[k]) {
        min[k] = al->x[i][k];
        *mini = i;
      }
      
      if(al->x[i][k] > max[k]) {
        max[k] = al->x[i][k];
        *maxi= i;
      }
      
        
    }
  }
  for(i=0;i<DIM;i++){
    box[i][i] = max[i] - min[i];
  }
}

/*===================================================*/


void do_fep_fit(t_atomlist *al, int f1, int f2, int f3,
                t_atomlist *al2, int g1, int g2, int g3)
{
  rvec ca1_;
  rvec ca2_;
  rvec diff_;
  rvec n1_;
  rvec n2_;
  rvec c1_;
  rvec c2_;
  int i,k;

  /* get vectors for fit atoms */

  copy_rvec(al->x[f1],ca1_);
  copy_rvec(al->x[f2],n1_);
  copy_rvec(al->x[f3],c1_);
  copy_rvec(al2->x[g1],ca2_);
  
  /* calculate vector from ca1 to ca2 */

  rvec_sub(ca1_,ca2_,diff_);
  
  /* correct atom position of al2 */

  for(k=0;k<al2->natoms;k++){
    rvec_add(al2->x[k],diff_,al2->x[k]);
  }

  /* get the vectors again */

  copy_rvec(al2->x[g1],ca2_);
  copy_rvec(al2->x[g2],n2_);
  copy_rvec(al2->x[g3],c2_);
  

  /*  make first rotation */
  /*    (all fit atoms in plane) */

  rvec ij;
  rvec ik;
  rvec is;
  rvec it;

  rvec_sub(ca1_,n1_,ij);
  rvec_sub(ca1_,c1_,ik);
  rvec_sub(ca1_,n2_,is);
  rvec_sub(ca1_,c2_,it);
  

  rvec ij_ik;
  rvec is_it;

  /* perpend. to the triangles */

  vectorProduct(ij,ik,ij_ik);
  vectorProduct(is,it,is_it);

  real absij_ik,absis_it;
  
  absij_ik = absVec(ij_ik);
  absis_it = absVec(is_it);

  real scp;
  real angle;

  scp = scalarProduct(ij_ik,is_it);
  angle = scp/(absij_ik*absis_it);
  angle=acos(angle);

  rvec rotvec;
  vectorProduct(ij_ik,is_it,rotvec);
  rvec rotsc;
  scaleVector(rotvec,rotsc);

  /* create the rotation matrices */

  matrix tm1;
  matrix tm2;
  
  createRotationMatrix1(rotsc,tm1);
  createRotationMatrix2(rotsc,tm2);

  /* rotate mutant around vector */

  for(i=0;i<al2->natoms;i++){
    rotateVectorAroundVector(al2->x[i],ca1_,tm1,tm2,-angle);
/*     al2=rotate_atom_around_bond(al2,tm1,tm2,-angle/180.0*PI,i,ca1_); */
  }

  /* get the vectors again */

  copy_rvec(al2->x[g1],ca2_);
  copy_rvec(al2->x[g2],n2_);
  copy_rvec(al2->x[g3],c2_);

  rvec_sub(ca1_,n2_,is);
  rvec_sub(ca1_,c2_,it);
  
  real absij,absis;

  scp=scalarProduct(ij,is);
  absij=absVec(ij);
  absis=absVec(is);
  angle=scp/(absij*absis);
  angle=acos(angle);
  scaleVector(ij_ik,rotsc);

  /* make new rotation matrices */

  createRotationMatrix1(rotsc,tm1);
  createRotationMatrix2(rotsc,tm2);

  /* determine direction to rotate */

  real d=determinant(ij,is,ij_ik);
  if(d>0) angle=-angle;

  /* make second rotation */

  for(i=0;i<al2->natoms;i++){
    rotateVectorAroundVector(al2->x[i],ca1_,tm1,tm2,angle);
  }

}


