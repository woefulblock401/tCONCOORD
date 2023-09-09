#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include "sysstuff.h" 
#include "typedefs.h"
#include "smalloc.h"
#include "copyrite.h"
#include "string2.h"
#include "macros.h"
#include "confio.h" 
#include "vec.h" 
#include "statutil.h" 
#include "futil.h"
#include "pdbio.h" 
#include "physics.h" 
#include "strdb.h" 
#include "gbutil.h" 
#include "index.h" 
#include "nsgrid.h"
#include "nrjac.h"
#include "mdatoms.h"
#include <tpxio.h>
#include "nrnb.h"
#include "ns.h"
#include "force.h" 
#include "constr.h"
#include "mdrun.h"
#include "disre.h"
#include "dihre.h"
#include "orires.h"
#include "do_fit.h"
#include "gmx_random.h"
#include "xtcio.h"
#include <structs.h>

#define MAKE_FULL_NEIGHBORLIST 0
#define UPDATE_NONBONDED 1



/*==========================================================*/

enum {SINGLE, DOUBLE, TRIPLE, CNSINGLE, CNOMEGA, CNDOUBLE, CNTRIPLE, NNSINGLE, NNDOUBLE, CODOUBLE, COSINGLE}; 

enum {BOND, B13, B14, NL};

enum {FREE, PHI, PSI, OMEGA, RING};

enum { rtOTH, rtPROTEIN, rtNUC, rtION, rtSOL};

enum {eNON, eCOV, eSUL, eBBHD, eBBHA, eSBH_BD, eSBH_BA, eSBH_SD, eSBH_SA,eSSHD,eSSHA, eNPR, ePHO, ePACK, eMET, eEX, eFO};

enum {ehNON, ehBBBB, ehBBSC, ehSCBB, ehSCSC };
  
enum {A, B, G, D, E, Z, H};


/*==========================================================*/

#define CCONTR   0.9
#define HCONTR   0.35
#define NCONTR   0.9
#define OCONTR   0.9
#define SCONTR   1.1
#define PCONTR   1.1
#define MCONTR   1.5
#define ICONTR   1.5
#define FCONTR   0.7
#define CLCONTR  1.0
#define BRCONTR  1.25

/*==========================================================*/


#define PI 3.141592653589793115997963468544

#define DIST(al,i,j) dist_ij(al->x[i],al->x[j])

#define DIST2(al,i,j) dist2(al->x[i],al->x[j])

#define DIHED(al,i,j,k,l) RAD2DEG*dihedral(al->x[i],al->x[j],al->x[k],al->x[l])


#define FPEqualTo(a,b)          (fabs((a)-(b)) < 1e-4)  
#define FPNotEqualTo(a,b)       (fabs((a)-(b)) > 1e-4)  
#define SetEulerAngles(ang1, ang2, ang3) \
        *Rx = (ang1);\
        *Ry = (ang2);\
        *Rz = (ang3)


/*==========================================================*/

#define IS_OMEGA(al,i,j) ((strcmp(al->name[i]," N  ") == 0 &&   \
                           strcmp(al->name[j]," C  ") == 0) ||\
                          (strcmp(al->name[i]," C  ") == 0 && \
                           strcmp(al->name[j]," N  ") == 0 )) 

#define IS_SSBOND(al,i,j) (strcmp(al->name[i]," SG ") == 0 &&\
                           strcmp(al->name[j]," SG ") == 0 && \
                           strcmp(al->resname[i],"CYS") == 0 && \
                           strcmp(al->resname[j],"CYS") == 0)

#define IS_PSI(al,i,j) (!strcmp(al->name[i]," CA ") && \
                           !strcmp(al->name[j]," C  ")) ||\
                          (!strcmp(al->name[i]," C  ") && \
                              !strcmp(al->name[j]," CA ")) 

#define IS_PHI(al,i,j) (!strcmp(al->name[i]," CA ") && \
                           !strcmp(al->name[j]," N  ")) ||\
                          (!strcmp(al->name[i]," N  ") && \
                              !strcmp(al->name[j]," CA ")) 


#define IS_HYDROPHOBIC(al,i) (!strcmp(al->symbol[i],"C") && \
                              (!strcmp(al->resname[i],"ILE") ||\
                               !strcmp(al->resname[i],"PHE") ||\
                               !strcmp(al->resname[i],"MET") ||\
                               !strcmp(al->resname[i],"LEU") ||\
                               !strcmp(al->resname[i],"TRP") ||\
                               !strcmp(al->resname[i],"VAL")) &&\
                               al->order[i] != 0)

#define IS_CHARGED_RESI(al,i) (!strcmp(al->resname[i],"ARG") || \
                               !strcmp(al->resname[i],"LYS") || \
                               !strcmp(al->resname[i],"GLU") || \
                               !strcmp(al->resname[i],"ASP"))

/*==========================================================*/


static char *pdbstr = "ATOM  %5d %-4s%1s%3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n";


static char *prot_res [] = {
  "ALA","ARG","ASN","ASP","ASPH","CYS","GLN","GLU","GLUH","GLY","HIS",
  "HISA","HISB","HISH","ILE","LEU","LYS","LYSH","MET","PHE","PRO","SER","THR","TRP",
  "TYR","VAL","HIE","HID"
};


static char *nucac_res [] = {
  "  A","  C","  G","  U","  T"
};

static char *ion_res [] = {
  " NA","  K"," MN"," FE"," CA"," MG",
  " CU"," ZN"," CL","SO4","SUL","PO4"
};

static char *sol_res [] = {
  "SOL","HOH","HO4"
};

/*==========================================================*/

static int  MAX_BOND =  10;
static int  MAX_B13  =  20;
static int  MAX_B14  =  40;
static int  MAX_NL   =  2000;
static int  MAX_CON  =  200;


/*==========================================================*/
/*                      FUNCTIONS                           */
/*==========================================================*/

/*==========================================================*/
/*                                                          */
/*              Structure initialisation                    */
/*                                                          */
/*==========================================================*/

extern t_atomlist *atomlist_init(void);

extern t_contab *contab_init(void);

extern t_bounds *bounds_init(int n);

extern t_boundtrack *boundtrack_init(void);

extern t_dihed *dihed_init(void);

extern t_resl *resl_init(void);

extern t_vdwcomb *vdwcomb_init(void);

extern t_types *types_init(void);

extern t_vdw *vdw_init(void);

extern t_bondlist *bl_init(void);

extern t_idxgroups *idx_init(void);

extern t_input *input_init(void);

extern t_excl *excl_init(void);

extern t_force *force_init(void);



/*==========================================================*/
/*                                                          */
/*              Structure reallocation                      */
/*                                                          */
/*==========================================================*/

extern t_atomlist *al_realloc(t_atomlist *al,int natoms);

extern t_contab *contab_realloc(t_contab *ct, int n);

extern t_bounds *bounds_realloc(t_bounds *b, int n);

extern t_boundtrack *boundtrack_realloc(t_boundtrack *bt,int n);

extern t_dihed *dihed_realloc(t_dihed *dh, int n);

extern t_resl *resl_realloc(t_resl *rl, int n);

extern t_idxgroups *idx_realloc(t_idxgroups *grp, int n);

extern t_vdwcomb *vdwcomb_realloc(t_vdwcomb *vdw, int n);

extern t_bondlist *bl_realloc(t_bondlist *bl, int n);

/*==========================================================*/
/*                                                          */
/*              Structure deallocation                      */
/*                                                          */
/*==========================================================*/

extern void free_al(t_atomlist *al);

extern void free_contab(t_contab *ct);

extern void free_groups(t_idxgroups *grp);

extern void free_vdwcomb(t_vdwcomb *vdw);

extern void free_types(t_types *tp);

extern void free_vdw(t_vdw *vdw);


/*==========================================================*/
/*                                                          */
/*              Error Messages                              */
/*                                                          */
/*==========================================================*/

extern void CNCwarn(FILE *log,char *warn);

extern void CNCerr(FILE *log,char *err);

extern void CNClog(FILE *log,char *string);

extern void fatal_error(char *error);

/*==========================================================*/
/*                                                          */
/*              atomlist.c                                  */
/*                                                          */
/*==========================================================*/

extern void copy_atom(t_atomlist *al1, int i, t_atomlist *al2, int k);

extern t_atoms *al2atoms(t_atoms *atoms,t_atomlist *al, rvec **x);

extern void rvec2al(t_atomlist *al,rvec *x);

extern void al2rvec(rvec *x,t_atomlist *al);

extern void get_symbol(FILE *log,t_atomlist *al);

extern void get_symbol2(FILE *log,t_atomlist *al);

extern void get_bconstr(t_atomlist *al);

extern void com(t_atomlist *al);

extern void calc_cent(t_atomlist *al, rvec com);

extern void com_frag(t_atomlist *al, int n,int *idx,rvec *new,rvec cm);

extern bool al_connected(t_atomlist *al, int i, int k);

extern int count_res(t_atomlist *al);

extern real get_vdw(t_atomlist *al, int i, int k, int flag);

extern void print_atom(FILE *fp,t_atomlist *al,int i);

extern void write_pdb(t_atomlist *al,char *filename);

extern void write_pdb_frame(FILE *fp,t_atomlist *al, int n);

extern void al2xtc(t_atomlist *al, int xtc, int frame);

extern void renumber_atoms(t_atomlist *al, int start);

extern void renumber_residues(t_atomlist *al, int start);

extern void rename_at(t_atomlist *al);

extern void get_order(t_atomlist *al);

extern void occ_to_one(t_atomlist *al);

extern t_atomlist *al_from_atoms(t_atoms *atoms,
                                 rvec *x);

extern void bonds_from_idef(t_atomlist *al, t_idef *idef);
extern t_atomlist *fuse_al(t_atomlist *al1, t_atomlist *al2);



/*==========================================================*/
/*                                                          */
/*              angle.c                                     */
/*                                                          */
/*==========================================================*/

extern void get_angle_vec(rvec x1,rvec x2, rvec x3, rvec shift);

extern int check_angle(rvec at3, rvec at1, rvec at2, t_bounds *b, 
                       int i, bool notouch1, bool notouch2);


/*==========================================================*/
/*                                                          */
/*              array.c                                     */
/*                                                          */
/*==========================================================*/

extern void array_sort(int *iarr, real *darr,int size);

extern void copy_iarr(int *dest, int *src, int size);

extern void i_sort(int *array, int size);

extern void switch_real(real *x, real *y);

extern void switch_int(int *x, int *y);

/*==========================================================*/
/*                                                          */
/*              bonded.c                                    */
/*                                                          */
/*==========================================================*/

extern real do_angle_force(t_atomlist *al, t_bounds *b, real weight);


extern real do_bound_force(t_atomlist *al, t_bounds *b, 
                           real bweight, real aweight,real dweight,
                           gmx_rng_t rng);

extern real do_bound_force2(t_atomlist *al, t_bounds *b, real weight);

/*==========================================================*/
/*                                                          */
/*              bondlib.c                                   */
/*                                                          */
/*==========================================================*/


extern t_bondlib *read_bonds_and_angles(int bonds);

extern bool get_bond(char *type1, char *type2, t_bondlib *bl, real *av, real *lb, real *ub);

extern bool get_angle(char *type1, char *type2, char *type3, t_bondlib *bl,
                      real *av, real *lb, real *ub, real *avang, real *sig);


/*==========================================================*/
/*                                                          */
/*              bondlist.c                                  */
/*                                                          */
/*==========================================================*/

extern void restrict_bond(t_bondlist *bl, int at1, int at2);

extern bool is_bond(t_atomlist *al, int at1, int at2);

extern void rotatable_bonds(t_atomlist *al,t_bondlist *bl);

extern void print_rotatable_bonds(FILE *log,t_atomlist *al, t_bondlist *bl);


/*==========================================================*/
/*                                                          */
/*              bounds.c                                    */
/*                                                          */
/*==========================================================*/


extern void input_error(int n, char *line);

extern int count_bounds(FILE *fp);

extern t_bounds *bounds_init(int n);

extern t_bounds *read_bounds(char *filename, bool bVerbose);

extern void copy_bound(t_bounds *b1, int i, t_bounds *b2, int k);

extern void check_bounds(FILE *log,t_bounds *b,t_atomlist *al, t_idxgroups *pln,
                         t_idxgroups *imp);

extern void damp_bounds(t_bounds *b, real damp);


/*==========================================================*/
/*                                                          */
/*              boundtrack.c                                */
/*                                                          */
/*==========================================================*/


extern void add_bound(t_boundtrack *bt, int i,int k);

extern void add_bound(t_boundtrack *bt, int i,int k);

extern bool is_bound(t_boundtrack *bt, int i, int k);

extern void fill_boundtrack(t_bounds *b, t_boundtrack *bt);


/*==========================================================*/
/*                                                          */
/*              buried.c                                    */
/*                                                          */
/*==========================================================*/
extern int find_atom(t_atomlist *al, rvec x, real cut2, int *exclude, int nex);

extern real ray_search(t_atomlist *al,rvec x, int *exclude, int nex);


/*==========================================================*/
/*                                                          */
/*              chiral.c                                    */
/*                                                          */
/*==========================================================*/

extern void get_impropers(t_atomlist *al, t_dihed *impr);

extern void print_improps(t_atomlist *al, t_dihed *impr);

extern void swap_sidechain(t_atomlist *al, t_resl *rl, int center, 
                           int *idx, int *tags);

extern int check_impr(t_atomlist *al, t_idxgroups *imp, int *imps);

extern int do_improp2(t_atomlist *al, t_dihed *impr, t_resl *rl,int *tags);

extern int do_improp(t_atomlist *al, t_idxgroups *imp, int *tags, 
                      gmx_rng_t rng);

extern int do_chiral(t_atomlist *al, t_idxgroups *imp, int *tags);

extern void mirror(t_atomlist *al);

extern real do_chiral_force(t_atomlist *al, t_idxgroups *imp, real weight);


/*==========================================================*/
/*                                                          */
/*              cnctop.c                                    */
/*                                                          */
/*==========================================================*/


extern void read_cnctop(char *filename, t_atomlist *al,t_idxgroups *acc, t_idxgroups *don,
                        t_idxgroups *phob, t_idxgroups *pl, t_idxgroups *imp);

  
extern void write_cnctop(char *filename,t_atomlist *al,t_idxgroups *don,
                         t_idxgroups *acc, t_idxgroups *phob, t_idxgroups *pln,
                         t_idxgroups *imp,t_input *inp, t_vdw *vdw, t_vdwcomb *vdwcomb,
                         int psize, atom_id *pos_id);

/*==========================================================*/
/*                                                          */
/*              contab.c                                    */
/*                                                          */
/*==========================================================*/


extern bool ct_connected(t_contab *ct, int i, int k, int flag1, int flag2);

extern void add2contab(t_contab *ct, int i, int k, int flag1, int flag2,bool check);

extern void add_bond(t_atomlist *al,int i, int k, t_contab *ct,int flag1, int flag2);

extern bool connected(t_contab *ct, int i, int k);

extern void print_contab(t_contab *rct, char *filename);

/*==========================================================*/
/*                                                          */
/*              dihed.c                                     */
/*                                                          */
/*==========================================================*/


extern void fill_dihed(t_dihed *dh, t_atomlist *al);

extern void get_dihedrals(t_atomlist *al, t_dihed *dihed);

/*==========================================================*/
/*                                                          */
/*              disco_util.c                                */
/*                                                          */
/*==========================================================*/

extern void copy_coords(rvec *source, rvec *target, int natoms);


extern int do_disco(FILE *log,t_atomlist *al, t_bounds *b,
                    t_boundtrack *bt,t_dihed *impr, t_idxgroups *imp,t_resl *rl, 
                    real bump_tol, int funcnr, int dfunc,
                    t_idxgroups *pl,gmx_rng_t rng, t_gridmap *gp,
/*                     t_topology *top, t_forcerec *fr, */
/*                     t_mdatoms *md, t_inputrec *ir, */
/*                     t_groups *grps, t_nrnb *nrnb, */
/*                     t_commrec *cr, t_block *cgs, t_nsborder *nsb, */
                    real rlong, int nl, int maxit, int impr_func,bool bVerbose,
                    bool bTarget,
                    t_atomlist *tal, int targ_n, atom_id *targ_index,
                    rvec *oldx, int cfreq,int npos, atom_id *posr,
                    real maxff, real nbfac, bool bDoRef, int nsteps, real gyr, real gyr_tol);

extern bool ligand_disco(t_atomlist *al, t_bounds *b, t_idxgroups *imp, t_idxgroups *pl,
                         int maxit, gmx_rng_t rng, FILE *log, t_boundtrack *bt, bool bVerbose);



extern int do_dock(FILE *log,t_atomlist *al, t_bounds *b, t_bounds *plb,
                   t_boundtrack *bt,t_dihed *impr, t_idxgroups *imp,t_resl *rl, 
                   real bump_tol, int funcnr, int dfunc,
                   t_idxgroups *pl,gmx_rng_t rng,
                   t_topology *top, t_forcerec *fr,
                   t_mdatoms *md, t_inputrec *ir,
                   t_groups *grps, t_nrnb *nrnb,
                   t_commrec *cr, t_block *cgs, t_nsborder *nsb,
                   real rlong, int nl, int maxit, int impr_func,bool bVerbose,
                   rvec *oldx, int cfreq, int npos, atom_id *posr);


/*==========================================================*/
/*                                                          */
/*              distance.c                                  */
/*                                                          */
/*==========================================================*/

extern void correct_dist(rvec x1, rvec x2, rvec diff, real d,
                         t_bounds *b, int i, int dfunc, gmx_rng_t rng,
                         bool notouch1, bool notouch2);



/*==========================================================*/
/*                                                          */
/*              dist_util.c                                 */
/*                                                          */
/*==========================================================*/

extern int write_bonds(FILE *fp, FILE *log,t_atomlist *al, t_boundtrack *bt, 
                       t_input *inp,t_excl *ex,t_bondlib *bl, bool bIgn);

extern int write_bonds_from_tpx(FILE *fp, FILE *log,t_atomlist *al, t_boundtrack *bt, 
                                t_input *inp,t_excl *ex,t_topology *top, bool bIgn);


extern int write_angles(FILE *fp, FILE *log,t_atomlist *al, t_boundtrack *bt, 
                        t_input *inp,t_excl *ex, t_bondlib *bl, bool bIgn);

extern int write_angles_from_tpx(FILE *fp, FILE *log,t_atomlist *al, 
                                 t_boundtrack *bt, t_input *inp,t_excl *ex, 
                                 t_topology *top, bool bIgn);

  
extern int write_hbond_angles(FILE *fp, FILE *log, t_atomlist *al, 
                              t_boundtrack *bt,t_excl *ex);

extern int write_dihedrals(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt, 
                           t_dihed *dihed,
                           t_input *inp,t_excl *ex);

extern int write_planar(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt, 
                        t_idxgroups *pln,t_input *inp,t_excl *ex);


extern int write_sidechain(FILE *fp, FILE *log, t_atomlist *al, t_resl *rl, 
                           t_boundtrack *bt, int min);


extern int write_hbonds(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt,
                         t_excl *ex);


extern int write_hydrophobics(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt, 
                              t_contab *rct, t_resl *rl, t_idxgroups *phb,
                              t_input *inp, gmx_rng_t rng, real def, t_excl *ex,
                              bool bnoH);

extern int write_close_pairs(FILE *fp, FILE *log, t_atomlist *al,t_boundtrack *bt,
                             t_input *inp,int min, real fixdist, t_excl *ex);

extern int write_packing_constraints(FILE *fp, FILE *log, t_atomlist *al,t_boundtrack *bt,
                                     t_excl *ex, real pack_limit, t_input *inp);


extern int write_network_restrictions(FILE *fp, FILE *log, t_atomlist *al, t_boundtrack *bt, 
                                      t_idxgroups *resgroups,
                                      t_resl *rl,int max, gmx_rng_t rng,
                                      t_input *inp, real def, int maxn,
                                      t_excl *ex, bool bnoH);

extern int write_forced(FILE *fp, FILE *log, t_atomlist *al, t_resl *rl, t_force *fo, 
                        t_input *inp, t_boundtrack *bt);


extern int long_range(FILE *fp, FILE *log, t_atomlist *al,t_boundtrack *bt,t_resl *rl, 
                      int max,gmx_rng_t rng, 
                      t_excl *ex, t_input *inp);

extern int write_neighbor_res(FILE *fp, FILE *log, t_atomlist *al,t_boundtrack *bt,
                              t_excl *ex, bool bIgn);

extern void print_constr_per_atom(FILE *fp, t_boundtrack *bt, t_atomlist *al);


/*==========================================================*/
/*                                                          */
/*              exfo.c                                      */
/*                                                          */
/*==========================================================*/


extern void read_exfo(char *filename, t_excl *ex, t_force *fo);

extern bool is_excluded(t_excl *ex, t_atomlist *al,int id1, int id2);

extern bool is_forced(t_force *fo, t_atomlist *al,int id1, int id2);


/*==========================================================*/
/*                                                          */
/*              force.c                                     */
/*                                                          */
/*==========================================================*/

extern void dump_forces(t_atomlist *al, char *filename);

extern void apply_f(t_atomlist *al, real weight, int *tags);
extern void apply_f2(t_atomlist *al, real weight, int *tags, real max);


extern void clear_forces(t_atomlist *al);

extern real sum_forces(t_atomlist *al);

extern void backbone_opt(t_atomlist *al, int *nt, int *tags);

/*==========================================================*/
/*                                                          */
/*              geometry.c                                  */
/*                                                          */
/*==========================================================*/

extern real dist2(rvec v1, rvec v2);
/* return squared distance */

extern real dist_ij(rvec v1, rvec v2);
/* return distance */
extern real torsion(rvec v1,rvec v2, rvec v3, rvec v4);
/* return the torsion angle */
extern real dihedral(const rvec xi,const rvec xj,const rvec xk,const rvec xl);
/* gromacs dihedral angle */

extern real angle_ij_ik(rvec v1,rvec v2, rvec v3);

extern real scalarProduct(rvec v1, rvec v2);
/* calculate scalar product */
extern real absVec(rvec v);

extern void vectorProduct(rvec v1, rvec v2, rvec v);

extern void vectorPlusVector(rvec v1,rvec v2,rvec v);

extern void vectorMinusVector(rvec v1,rvec v2,rvec v);

extern real determinant(rvec a,rvec b, rvec c);

extern void scaleVector(rvec v_, rvec vsc_);

extern void createRotationMatrix1(rvec a_,matrix m);

extern void createRotationMatrix2(rvec a_, matrix m);

extern void MatrixVectorMultiply(matrix m,rvec x_, rvec v_);

extern void rotateVectorAroundVector(rvec v1_,rvec v2_,matrix tm1,
                              matrix tm2, real phi);

extern void rotate_atom(rvec v1_, rvec v2_, rvec v3_, real phi);

extern void princ_comp(int n, rvec x[],
                    matrix trans,rvec d);

extern void rotate_atoms(int n,rvec x[],matrix trans);

extern void max_coords2(t_atomlist *al, matrix box, int *maxi, int *mini);

extern void max_coords(t_atomlist *al, matrix box);

extern real calc_rmsd(rvec *x1, rvec *x2, int size);

extern int findEulerAngles(matrix3x3 R, double *Rx, double *Ry, double *Rz);

extern void do_fep_fit(t_atomlist *al, int f1, int f2, int f3,
                       t_atomlist *al2, int g1, int g2, int g3);



/*==========================================================*/
/*                                                          */
/*              gridmap.c                                    */
/*                                                          */
/*==========================================================*/



extern t_gridmap *read_map(char *filename, FILE *log,bool bVerbose);

extern void write_dx_map(char *filename, t_gridmap *map);

extern t_gridmap *gridmap_init(void);

extern t_gridmap *reset_grid(t_gridmap *map);

/*==========================================================*/
/*                                                          */
/*              groups.c                                    */
/*                                                          */
/*==========================================================*/

extern void add_to_group(t_idxgroups *grp, int ir, int id);
extern void add_to_group_nosort(t_idxgroups *grp, int ir, int id);


extern t_idxgroups *del_group(t_idxgroups *grp, int i);

extern void sort_group_by_order(t_atomlist *al, t_idxgroups *idx);

extern void sort_group_by_posres(t_atomlist *al, t_idxgroups *idx);

extern bool is_same_group(t_idxgroups *grp, int i, int k);

extern bool is_redundant_group(t_idxgroups *grp, int i, int k);

extern void print_groups(t_atomlist *al,t_idxgroups *grp);

extern void print_group(t_atomlist *al,int *idx, int n);

extern void copy_group(t_idxgroups *grp1,int idx1,t_idxgroups *grp2,int idx2);

extern t_idxgroups *remove_redundant_groups(t_idxgroups *idx);

extern bool is_planar_group(t_atomlist *al,int n,int *atoms, 
                            const real tol, real *plan);

extern t_idxgroups *merge_groups(t_idxgroups *grp, int min, bool bVerbose);

extern void add_to_group_blind(t_idxgroups *grp, int ir, int id);

extern int is_in_group(t_idxgroups *grp, int idx, int item);

extern void sort_groups(t_idxgroups *grp);

extern void write_group_script(char *filename, char *sel, t_idxgroups *grp);


/*==========================================================*/
/*                                                          */
/*              hbonds.c                                    */
/*                                                          */
/*==========================================================*/

extern t_hbond *hbond_init(void);

extern t_hbond *copy_hbond(t_hbond *hb);
extern void free_hbonds(t_atomlist *al);

extern void add_hbond(t_atomlist *al, int don, int acc, int type, real pr_rad);

extern void get_hbonds(t_atomlist *al, real dmax, real amin, real pr_rad);

extern real hbond_energy(real d);

extern real hbond_protection(t_atomlist *al, int don, int acc, real rad);

extern void print_hbonds_to_file(FILE *fp, t_atomlist *al);

extern real energy_from_histogram(t_histogram *h, real x);


/*==========================================================*/
/*                                                          */
/*              histogram.c                                 */
/*                                                          */
/*==========================================================*/


extern t_histogram *new_histogram(real lb, real ub, real delta);

extern void free_histogram(t_histogram *h);

extern void add_value(t_histogram *h, real val, bool flag, real weight);

extern void print_histogram(FILE *fp,t_histogram *h);

extern void norm_histogram(t_histogram *h);

extern real integrate(t_histogram *h);

extern t_histogram *read_histogram(FILE *fp);


/*==========================================================*/
/*                                                          */
/*             input.c                                      */
/*                                                          */
/*==========================================================*/

extern void print_input_file(FILE *fp,t_input *inp);
extern void print_input_log(FILE *fp,t_input *inp);

extern t_input *read_cnc_input(FILE *log,char *filename);


/*==========================================================*/
/*                                                          */
/*             lib_util.c                                   */
/*                                                          */
/*==========================================================*/

extern FILE *cnclib(char *filename);

extern void get_cnc_solv(t_atomlist *al);

extern void get_types(FILE *log,t_atomlist *al, t_types *tp);

extern void get_vdw_radii(FILE *log,t_atomlist *al, t_vdw *vdw, 
                          t_vdwcomb *vdwcomb, t_vdwcomb *vdw14comb);

extern t_area *read_area(char *filename);


/*==========================================================*/
/*                                                          */
/*             namegroups.c                                 */
/*                                                          */
/*==========================================================*/


extern t_namegroups *read_namegroups(FILE *fp);


/*==========================================================*/
/*                                                          */
/*             nbsearch.c                                   */
/*                                                          */
/*==========================================================*/


extern void nb_search(t_atomlist *al,
                      t_contab *ct,real rlong, bool bVerbose);

extern void dump_nl(FILE *fp,t_atomlist *al,real max);

extern void update_neighborlist(FILE *log,t_atomlist *al, rvec *x, t_topology *top,
                                t_mdatoms *md, t_forcerec *fr,
                                t_inputrec *ir,
                                t_groups *grps, t_nrnb *nrnb,
                                t_commrec *cr, t_block *cgs,
                                t_nsborder *nsb, real rlong, bool bSimple,
                                int nat, atom_id *ids);

extern void update_nl(t_atomlist *al, t_nblist *nlist);

extern void fill_nl(t_atomlist *al,t_nblist *nlist,t_contab *ct);

extern t_gridmap *nb(t_atomlist *al, t_gridmap *gp, real cutoff);
extern int gridp(real x, real origin, real inv_spacing, int max);

extern void fill_neighborlist(FILE *log, t_atomlist *al, t_gridmap *gp, real cutoff, bool bUpdate);

extern void free_gridmap(t_gridmap *gp);


/*==========================================================*/
/*                                                          */
/*             nonbonded.c                                  */
/*                                                          */
/*==========================================================*/


extern real do_nb_force(t_atomlist *al, t_boundtrack *bt, 
                        real tol, real weight);


extern int check_bumps(t_atomlist *al, t_boundtrack *bt,
                t_vdwcomb *vdwcomb, real tol, gmx_rng_t rng, 
                       int funcnr, int *tags, bool bGroup, int grpnr);

extern int count_bumps(t_atomlist *al, t_boundtrack *bt, real tol, real *bsum, real max, bool bVerbose);

extern void bumps_check(FILE *of,t_atomlist *al, t_boundtrack *bt, 
                        real tol);

/*==========================================================*/
/*                                                          */
/*             packing.c                                    */
/*                                                          */
/*==========================================================*/
extern real atom_packing(t_atomlist *al, int i);
extern real packing_score(t_atomlist *al);
extern real packing_score_subset(t_atomlist *al, atom_id *atoms, int size);
  
extern int count_contacts(t_atomlist *al, atom_id *atoms, int size, real cutoff);



/*==========================================================*/
/*                                                          */
/*             pdbqt.c                                     */
/*                                                          */
/*==========================================================*/

extern t_atomlist *read_pdbqt(char *filename, FILE *log, bool bVerbose);


/*==========================================================*/
/*                                                          */
/*             planar.c                                     */
/*                                                          */
/*==========================================================*/

extern bool check_planar(t_atomlist *al, int n, int *idx, real tol);

extern bool count_planar(t_atomlist *al, t_idxgroups *pl, real *sum,
                         int *count, real max, real *pworst);

extern real do_planar_force(t_atomlist *al, int n, int *idx, real tol, real weight);

extern int get_planar_groups(FILE *log, t_atomlist *al, t_resl *rl, 
                              t_idxgroups *pln, bool bIgn);

/*==========================================================*/
/*                                                          */
/*             rama.c                                     */
/*                                                          */
/*==========================================================*/

extern t_gridmap *read_rama(void);
extern real val_from_2d_table(t_gridmap *gp, real v1, real v2);

/*==========================================================*/
/*                                                          */
/*             random.c                                     */
/*                                                          */
/*==========================================================*/

extern int random_int(gmx_rng_t rng, int max);

extern int random_int_in_range(gmx_rng_t rng, int min, int max);

extern void rand_array_int(gmx_rng_t rng, int n, int *array);

extern void simple_array_int(int *array, int n);

extern void random_al_simple(t_atomlist *al, gmx_rng_t rng);

extern void random_al(t_atomlist *al,gmx_rng_t rng, bool bTarget, rvec tco);

extern void perturb_al(t_atomlist *al,gmx_rng_t rng);

extern void random_al_idx(t_atomlist *al,int size, atom_id *index, gmx_rng_t rng);

extern void random_al_idx2(t_atomlist *al,int size, atom_id *index, 
                           gmx_rng_t rng, rvec cent, rvec bsize);


extern void randomise_group(t_atomlist *al, gmx_rng_t rng, int grpnr);


/*==========================================================*/
/*                                                          */
/*             resl.c                                       */
/*                                                          */
/*==========================================================*/


extern void fill_resl(t_resl *rl,t_atomlist *al);
extern int get_atom_idx(t_resl *rl, int residx, t_atomlist *al, char *name);



/*==========================================================*/
/*                                                          */
/*             rosetta.c                                       */
/*                                                          */
/*==========================================================*/

extern void get_rosetta_types(t_atomlist *al);
extern real lj_lattice_energy(t_atomlist *al,rvec x,real cutoff, 
                              real eps, real r, int maxorder);

extern real rs_solv_energy(t_atomlist *al);
extern real rs_hbond_energy(t_atomlist *al);
extern real rs_lj_energy(t_atomlist *al);
extern real bb_hbond_energy(t_atomlist *al);

/*==========================================================*/
/*                                                          */
/*             rotation.c                                   */
/*                                                          */
/*==========================================================*/

extern t_rotation *rotation_init(void);
extern t_rotation *rotation_realloc(t_rotation *rot, int n);
extern t_rotdata *rotdata_init(void);
extern t_rotdata *rotdata_realloc(t_rotdata *rt, int n);
extern t_rotdata *read_rotation_lib(void);
extern void write_rotations(FILE *fp,t_atomlist *al, t_resl *rl, t_rotdata *rt);
extern t_idxgroups *read_rotations(char *filename, bool bVerbose);
extern void rotate_group(t_atomlist *al, t_idxgroups *rot, int idx, real phi);
extern void do_random_rotations(t_atomlist *al, t_idxgroups *rot, int n, gmx_rng_t rng);
extern void make_bb_rotations(FILE *fp, t_atomlist *al, t_resl *rl, int nrot, bool bVerbose);

/*==========================================================*/
/*                                                          */
/*             sasa.c                                  */
/*                                                          */
/*==========================================================*/

extern void get_sasa_types(t_atomlist *al);
extern real sasa(t_atomlist *al);
extern real sasa_energy(t_atomlist *al);


/*==========================================================*/
/*                                                          */
/*             sidechain.c                                  */
/*                                                          */
/*==========================================================*/

extern t_scdata *read_scdata(void); 


/*==========================================================*/
/*                                                          */
/*             string_util.c                                */
/*                                                          */
/*==========================================================*/

extern void cut_string(char *in, char c);

extern void slice_string(char *dest, char *src, int beg, int end);

extern void switch_string(char *s1,char *s2);

extern bool find_key_word(FILE *fp,char *key);

extern void strip(char *out,char *in);

extern void substring(char *out, char *in, char begin, char end);

extern void insert_string(char *string, char *insert, int pos, char *flag);

extern void extend_name(char *name);

extern void progress(FILE *fp,char *string, int now, int full);

extern void cnc_copyright(char *prgname);

extern int str_nelem(char *str,int maxptr,char *ptr[]);

extern void papers_log(FILE *log);


/*==========================================================*/
/*                                                          */
/*             structure.c                                  */
/*                                                          */
/*==========================================================*/


extern bool is_polar_H(t_atomlist *al, int i);

extern bool IsProtein(t_atomlist *al, int k);

extern bool IsNucac(t_atomlist *al, int k);

extern bool IsIon(t_atomlist *al, int k);

extern bool IsSol(t_atomlist *al, int k);

extern void get_ring(t_atomlist *al, int i, t_idxgroups *tr, t_bondlist *bl);

extern void hbonds(t_atomlist *al, t_idxgroups *don, t_idxgroups *acc);

extern void fake_charge(t_atomlist *al);

extern void hydrophobics(t_atomlist *al, t_idxgroups *phob, FILE *fp);

extern void do_cov(t_atomlist *al, t_bondlist *bl, t_contab *rct, t_excl *ex);

extern void do_hbonds(t_atomlist *al, t_contab *rct,
                      FILE *fp,bool bCheck[], 
                      real hbdist,real hbangle, 
                      rvec minHP, real *mindist,
                      int *bbhb, int *bbhb_rej,
                      int *scbbhb, int *scbbhb_rej,
                      int *scschb, int *scschb_rej);

extern void do_hydrophobics(t_atomlist *al,t_contab *rct,real min,
                            t_idxgroups *phb);
extern void print_hydrophobics_to_file(FILE *fp, t_atomlist *al, 
                                       t_idxgroups *phb);


extern void get_hybrid(t_atomlist *al);

extern void get_hybrid2(t_atomlist *al, t_types *tp);

extern real radius_of_gyration(t_atomlist *al, rvec gvec);

extern void sequence(t_resl *rl, char *seq);

extern void print_sequence_log(FILE *fp, char *seq, int n);

/*==========================================================*/
/*                                                          */
/*             types.c                                      */
/*                                                          */
/*==========================================================*/

extern t_types *read_atom_types(int num);


/*==========================================================*/
/*                                                          */
/*             vdw.c                                        */
/*                                                          */
/*==========================================================*/


extern t_vdw *read_vdw_radii(int num);


/*==========================================================*/
/*                                                          */
/*             vdwcomb.c                                    */
/*                                                          */
/*==========================================================*/

extern t_vdwcomb *read_vdwcomb(int num, bool flag);

extern t_vdwcomb *read_vdwcomb2(FILE *fp, bool flag);



















