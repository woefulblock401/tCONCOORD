/*==========================================================*/
/* CONCOORD STRUCTURES AND TYPE DEFINITIONS                 */
/*==========================================================*/

#include<typedefs.h>

typedef char shstr[10];
typedef char longstr[256];
typedef double matrix3x3 [3][3];

typedef struct t_atomlist t_atomlist;
typedef struct t_pdbrecord t_pdbrecord;
typedef struct t_contab t_contab;
typedef struct t_dihed t_dihed;
typedef struct t_resl t_resl;
typedef struct t_vdwcomb t_vdwcomb;
typedef struct t_vdw t_vdw;
typedef struct t_list t_list;
typedef struct t_types t_types;
typedef struct t_namegroups t_namegroups;
typedef struct t_idxgroups t_idxgroups;
typedef struct t_bondlist t_bondlist;
typedef struct t_excl t_excl;
typedef struct t_force t_force;
typedef struct t_bounds t_bounds;
typedef struct t_boundtrack t_boundtrack;
typedef struct t_bondlib t_bondlib;
typedef struct t_input t_input;
typedef struct t_area t_area;
typedef struct t_histogram t_histogram;
typedef struct t_hbond t_hbond;
typedef struct t_scdata t_scdata;
typedef struct t_gridmap t_gridmap;
typedef struct t_rotation t_rotation;
typedef struct t_rotdata t_rotdata;


/*==========================================================*/
struct t_atomlist 

/* list of atoms containing all atom properties */

{
  int natoms;              /* number of atoms */
  int ngroups;             /* number of groups */
  int *id;                 /* array of atom id's */
  int *resid;              /* array of residue id's */
  shstr *name;             /* array of atom names */
  shstr *resname;          /* array of residue names */
  shstr *chain;            /* array of chain names */
  shstr *type;             /* array of atom types */
  shstr *altloc;           /* array of alternate locations */
  rvec *x;                 /* array of coordinate vectors  */
  rvec *f;                 /* array of fake force vectors  */
  real *occ;               /* array of occupancies  */
  real *bfac;              /* array of b-factors    */
  shstr *symbol;           /* array of element symbols */
  shstr *hyb;              /* array of hybridisation strings */
  int *order;              /* array of atom orders */
  int **bonds;             /* array of bond arrays */
  int *nbonds;             /* number of bonds for each atom */
  int **b13;               /* array of 1-3 neighbors */
  int *nb13;               /* number of 1-3 neighbors for each atom */
  int **b14;               /* array of 1-4 neighbors */
  int *nb14;               /* number of 1-4 neighbors for each atom */
  int **nb;                /* array of neighbors */
  int *nnb;                /* number of neighbors for each atom */
  int *memcount;           /* array of memory counters */ 
  int *restype;            /* array of residue types */
  int *ngrps;              /* number of groups for each atom */
  int **grpnr;             /* array of group numbers */
  bool *isdon;             /* hbond donor flag; TRUE if atom is a hbond donor */
  bool *isacc;             /* hbond acceptor flag; TRUE if atom is hbond acceptor */
  bool *ishphob;           /* hydrophobicity flag; TRUE if atom is hydrophobic */
  bool *isposres;          /* position restraints flag; TRUE if atom position is constant */
  bool *isflex;            /* flexibility flag; TRUE if atom is flexible */
  real *bcontr;            /* constribution to a covalent bond */
  real *vdw;               /* van der Waals radius */
  real *vdw14;             /* van der Waals radius for 1-4 pairs */
  real **vdwtab;           /* vdw look up table */
  int nvdw;               /* dimension of vdwtab */
  real **vdw14tab;         /* vdw14 look up table */
  int *ptype;              /* shortcut for atomtype */
  real *m;                 /* mass */
  real *q;                 /* charge */
  real *cnc_solv;          /* CONCOORD solvation parameter */
  int *rs_type;            /* ROSETTA atom types */
  real *rs_rad;            /* ROSETTA atom radius */
  real *rs_G;              /* ROSETTA delta Gfree */
  real *rs_Gref;           /* ROSETTA delta Gref */
  real *rs_V;              /* ROSETTA atom volume */
  real *rs_eps;            /* ROSETTA epsilon */
  real *rs_lambda;         /* ROSETTA lambda */
  real **rs_comb;          /* ROSETTA radius combinations */
  int nrs_comb;            /* number of ROSETTA combinations */
  t_hbond **hbonds;
  int nhbonds;
  ivec *cell;
  int *sasa_type;
  real *sasa_Ri;
  real *sasa_pi;
  real *sasa_sigi;
};

/*==========================================================*/
struct t_pdbrecord
/* structure for reading pdb files */
{
  int n;              /* number of atomlists */
  t_atomlist **al;    /* array of atomlist */
};
/*==========================================================*/

struct t_contab
/* connection table structure */
{
  int n;         /* number of entries */
  int *ncon;     /* number of connections for each entry */
  int **con;     /* array of connection id's */
  int **type;    /* array of connection types */
};
/*==========================================================*/
struct t_dihed
/* structure for storing dihedrals and impropers*/
{
  int n;                     /* number of dihedrals */
  int *center;               /* array of center atoms (for impropers ) */
  int *at1,*at2,*at3,*at4;   /* array's of atom id's */
  shstr *tp1,*tp2,*tp3,*tp4; /* array's of atom types */
  real *phi;                 /* array of angles */
  int *flag;                 /* array of flags (OMEGA,PHI,PSI...) */
  bool *restr;               /* array of flags; TRUE if dihedral is rotatable */
};

/*==========================================================*/
struct t_resl
/* residue list structure */
{
  int nres;           /* number of residues */
  int *id;            /* array of residue id's */
  shstr *resname ;    /* array of residue names */
  int *j0,*j1;        /* array of atom id's to 
                         find the atoms in the 
                         corresponding atomlist 
                      */
  bool *phi_restr;    /* TRUE if phi angle is restricted */
  bool *psi_restr;    /* TRUE if psi angle is restricted */
  bool *sc_restr;     /* TRUE if the side-chain is restricted */
  int *ngrps;         /* number of groups */
  int **grps;         /* array of group id's */
};

/*==========================================================*/
struct t_vdwcomb
/* structure for storing van-der-Waals combinations */
{
  int n;          /* number of combinations */
  shstr *type1;   /* array of atom types */
  shstr *type2;   /* array of atom types */
  real *vdw;      /* array of distances  */
  int *rs_type1;
  int *rs_type2;
  real *rs_comb;
};
/*==========================================================*/
struct t_vdw
/* structure for assignment of van-der-Waals radii */
{
  int n;           /* number of entries */
  shstr *type;     /* array of atom types */
  real *vdw;       /* array of vdw radii  */
  real *vdw14;     /* array of vdw 1-4 radii */
  int *rs_type;
  real *rs_rad;
  real *rs_G;
  real *rs_Gref;
  real *rs_V;
  real *rs_lambda;
  real *rs_eps;
};

/*==========================================================*/
struct t_list
/* structure for file parsing */
{
  int n;            /* number of lines */
  longstr *lines;   /* array of lines */
};

/*==========================================================*/
struct t_types
{
  int n;
  shstr *resname;
  shstr *name, *type, *hyb;
  int *rs_type;

};

/*==========================================================*/
struct t_namegroups
{
  int n;
  int *natoms;
  shstr *resname;
  shstr **atomnames;
  real *val;
};
/*==========================================================*/
struct t_idxgroups
{
  int n;          /* number of groups */
  int *natoms;    /* size of individual groups */
  int **atoms;    /* array of atom indices */
  real *val;      /* value for each group  */
  bool *flag;     /* flag for each group */
};

/*=============================================================*/
struct t_bondlist

  /* bond structure. 
     contains information about bonds
   */
{
  int n;             /* number of bonds */
  int *at1;          /* atom id */
  int *at2;          /* atom id */
  real *blen;        /* bond length */
  int *type;         /* bond type */
  bool *restricted;  /* TRUE if bond is not rotatable */
};

/*=============================================================*/

struct t_excl
{
  int n;
  int *id1;
  int *id2;
};

/*==========================================================*/

struct t_force
{
  int n;
  int *id1;
  int *id2;
  real *lb;
  real *ub;
};

/*==========================================================*/

struct t_bounds
{
  int n;
  int *at1, *at2, *at3, *at4;
  real *av, *lb, *ub, *ang, *dih, *sig;
  bool *isang, *isdih, *isbond;
  bool *pdla, *pald, *phplhp;
  bool *ac_ac, *don_ac, *don_don;
  bool *bCheck;
};

/*==========================================================*/

struct t_boundtrack
{
  int natoms;
  int *id;
  int *n;
  int **con;
};

/*==========================================================*/

struct t_bondlib
{
  int n;
  shstr *type1;
  shstr *type2;
  shstr *type3;
  real *av;
  real *lb;
  real *ub;
  real *ang;
  real *sig;
  bool *is_bond;
  bool *is_ang;
};

/*==========================================================*/

struct t_input
{
  real bond_tol;
  real angle_tol;
  real angle_sig;
  real bump_tol;
  real bump_14tol;
  real plan_tol;
  real restr_dihed_tol;
  real non_restr_dihed_tol;
  bool use_hbond_angles;
  bool use_hbonds;
  bool use_sidechain;
  bool use_hydrophobics;
  bool use_network;
  bool use_close_pairs;
  bool use_packing;
  real pack_limit;
  real close_pairs_fix_dist;
  bool use_neighbor_ca;
  bool use_long_range;
  real lr_tol[2];
  real network_tol[2];
  int min;
  int max;
  real hbond_max_dist;
  real hbond_min_angle;
  bool use_hbond_protection[3];
/*   rvec radius_of_gyration; */
/*   rvec radius_of_gyration_tol; */
  rvec solvation_rad;
  rvec solvation_max;
  real hphob_dist;
};

/*==========================================================*/

struct t_area
{
  real xmin,xmax;
  real ymin,ymax;
  real zmin,zmax;
};

/*==========================================================*/

struct t_histogram
{
  int n;
  real *y;
  real *x;
  real delta;
  real inv_delta;
  real lb;
  real ub;
};

/*==========================================================*/

struct t_hbond
{
  int don;
  int acc;
  int don_b;
  int acc_b;
  int type;
  real d_a_dist;
  real d_ab_dist;
  real db_a_dist;
  real db_ab_dist;
  real angle;
  real angle2;
  real dihed;
  real energy;
  real prot;
  bool constr;
  real packing;
};

/*==========================================================*/

struct t_scdata
{
  int n;
  shstr *resname;
  shstr *name1;
  shstr *name2;
  real *lb;
  real *ub;
};

/*==========================================================*/
struct t_gridmap
{
  char name[STRLEN];
  ivec npts;
  ivec n;
  rvec center;
  rvec origin;
  int nelem;
  real spacing;
  real inv_spacing;
  real *values;
  int **cell;
  int *natom;
  char datafile[STRLEN];
  char molecule[STRLEN];
  char paramfile[STRLEN];
  real precision;
};

/*==========================================================*/
struct t_rotation
{
  int n;
  shstr resname;
  longstr *names;
};

/*==========================================================*/
struct t_rotdata 
{
  int n;
  shstr *resname;
  t_rotation **rots;
};
