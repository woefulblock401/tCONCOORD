#include <tconcoord.h>

real atom_packing(t_atomlist *al, int i)
{
  
  int k;
  real vdwr,vol,bumps,contacts,natoms,di,xx;
  bumps = contacts = 0.;
  natoms = (real) al->natoms;
  real ref=0.120769;
  for(k=0;k<al->nnb[i];k++){
    if(i<al->nb[i][k]){
      vdwr = get_vdw(al,i,al->nb[i][k],FALSE);
      di = DIST(al,i,al->nb[i][k]);
      xx = di/vdwr;
      if(xx < 1.0){
        vol = pow(vdwr,3)-pow(vdwr*0.8,3);
        bumps +=1./vol;
      }
      else if (xx>=1.0 && xx <=1.226){
        vol = pow(vdwr*1.226,3)-pow(vdwr,3);
        contacts+=1./vol;
      }
    }
  }
  return ((contacts-bumps))/ref;
}

real packing_score(t_atomlist *al)
{
  int i,k;
  real vdwr,vol,bumps,contacts,natoms,di,xx;
  bumps = contacts = 0.;
  natoms = (real) al->natoms;
  real ref=0.120769;
  
  for(i=0;i<al->natoms;i++){
    for(k=0;k<al->nnb[i];k++){
      if(i<al->nb[i][k]){
        vdwr = get_vdw(al,i,al->nb[i][k],FALSE);
        di = DIST(al,i,al->nb[i][k]);
        xx = di/vdwr;
        if(xx < 1.0){
          vol = pow(vdwr,3)-pow(vdwr*0.8,3);
          bumps +=1./vol;
        }
        else if (xx>=1.0 && xx <=1.226){
          vol = pow(vdwr*1.226,3)-pow(vdwr,3);
          contacts+=1./vol;
        }
      }
    }
  }
  return ((contacts-bumps)/natoms)/ref;
}

real packing_score_subset(t_atomlist *al, atom_id *atoms, int size)
{
  int i,k;
  int idx;
  
  real vdwr,vol,bumps,contacts,natoms,di,xx;
  bumps = contacts = 0.;
  natoms = (real) size;
  real ref=0.120769;
  
  for(i=0;i<size;i++){
    idx = atoms[i];
/*     printf("name = %s\n",al->name[idx]); */
    
    for(k=0;k<al->nnb[idx];k++){
      
    
      vdwr = get_vdw(al,idx,al->nb[idx][k],FALSE);
      di = DIST(al,idx,al->nb[idx][k]);
      xx = di/vdwr;
/*       printf("\tname = %s xx = %g\n",al->name[al->nb[idx][k]],xx); */
      if(xx < 1.0){
        vol = pow(vdwr,3)-pow(vdwr*0.8,3);
        bumps +=1./vol;
      }
      else if (xx>=1.0 && xx <=1.226){
        vol = pow(vdwr*1.226,3)-pow(vdwr,3);
        contacts+=1./vol;
      }
    }
  }
  return ((contacts-bumps)/natoms)/ref;
}
  
int count_contacts(t_atomlist *al, atom_id *atoms, int size, real cutoff)
{
  int i,k,idx;
  real di;
  int count = 0;
  for(i=0;i<size;i++){
    idx = atoms[i];
    for(k=0;k<al->nnb[idx];k++){
      di = DIST(al,idx,al->nb[idx][k]);
      if(di<cutoff){
        count++;
      }
    }
  }
  return count;
}

          

  
