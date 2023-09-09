#include <tconcoord.h>

/*==========================================================*/
t_input *input_init(void)
{
  t_input *inp=NULL;
  snew(inp,1);
  inp->bond_tol = 0.04;
  inp->angle_tol = 0.2;
  inp->angle_sig = 9.;
  inp->bump_tol = 0.1;
  inp->bump_14tol = 0.08;
  inp->plan_tol = 0.03;
  inp->restr_dihed_tol = 0.1;
  inp->non_restr_dihed_tol = 0.6;
  inp->use_hbond_angles = TRUE;
  inp->use_hbonds = TRUE;
  inp->use_sidechain = FALSE;
  inp->use_packing = FALSE;
  inp->pack_limit = .8;
  inp->use_hydrophobics = TRUE;
  inp->use_network = TRUE;
  inp->use_close_pairs = FALSE;
  inp->close_pairs_fix_dist = 5.;
  inp->use_long_range = TRUE;
  inp->min = 30;
  inp->max = 60;
  inp->lr_tol[0] = 0.5;
  inp->lr_tol[1] = 3.;
  inp->network_tol[0] = 0.9;
  inp->network_tol[1] = 1.1;
  inp->hbond_max_dist = 2.6;
  inp->hbond_min_angle = 120.;
  inp->use_hbond_protection[0] =
    inp->use_hbond_protection[1] =
    inp->use_hbond_protection[2] = TRUE;

/*   inp->radius_of_gyration[XX] =  */
/*   inp->radius_of_gyration[YY] =  */
/*     inp->radius_of_gyration[ZZ] = 0.; */
/*   inp->radius_of_gyration_tol[XX] =  */
/*   inp->radius_of_gyration_tol[YY] =  */
/*     inp->radius_of_gyration_tol[ZZ] = 0.; */
  
  
  inp->solvation_rad[0] =
    inp->solvation_rad[1] =
    inp->solvation_rad[2] = 6.;
  inp->solvation_max[0] = 
    inp->solvation_max[1] = 
    inp->solvation_max[2] = 2.2; 
  inp->hphob_dist = 6.;
  return inp;
}

/*==========================================================*/
void print_input_log(FILE *fp,t_input *inp)
{
  shstr dummy[3];
  int i;
  fprintf(fp,"tCNC__log_> bond-tol               = %4.2f\n",inp->bond_tol);
  fprintf(fp,"tCNC__log_> angle-tol              = %4.2f\n",inp->angle_tol);
  fprintf(fp,"tCNC__log_> angle-sig              = %4.2f\n",inp->angle_sig);
  fprintf(fp,"tCNC__log_> plan-tol               = %4.2f\n",inp->plan_tol);
  fprintf(fp,"tCNC__log_> bump-tol               = %4.2f\n",inp->bump_tol);
  fprintf(fp,"tCNC__log_> bump-14tol             = %4.2f\n",inp->bump_14tol);
  fprintf(fp,"tCNC__log_> restr-dihed-tol        = %4.2f\n",inp->restr_dihed_tol);
  fprintf(fp,"tCNC__log_> non-restr-dihed-tol    = %4.2f\n",inp->non_restr_dihed_tol);
  if(inp->use_hbonds) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"tCNC__log_> use-hbonds             = %-4s\n",dummy[0]);
  if(inp->use_hbond_angles) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"tCNC__log_> use-hbond-angles       = %-4s\n",dummy[0]);
  if(inp->use_sidechain) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"tCNC__log_> use-sidechain          = %-4s\n",dummy[0]);
  if(inp->use_hydrophobics) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"tCNC__log_> use-hydrophobics       = %-4s\n",dummy[0]);
  if(inp->use_network) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"tCNC__log_> use-network            = %-4s\n",dummy[0]);
  fprintf(fp,"tCNC__log_> network-tol            = %4.2f %4.2f\n",inp->network_tol[0],inp->network_tol[1]);
  if(inp->use_packing) strcpy(dummy[0],"yes"); 
  else strcpy(dummy[0],"no"); 
  fprintf(fp,"tCNC__log_> use-packing            = %-4s\n",dummy[0]);
  fprintf(fp,"tCNC__log_> pack-limit             = %4.2f\n",inp->pack_limit); 

  if(inp->use_close_pairs) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"tCNC__log_> use-close-pairs        = %-4s\n",dummy[0]);
  fprintf(fp,"tCNC__log_> close-pairs-fix-dist   = %4.2f\n",inp->close_pairs_fix_dist);
  if(inp->use_neighbor_ca) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"tCNC__log_> use-neighbor-ca        = %-4s\n",dummy[0]);
  if(inp->use_long_range) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"tCNC__log_> use-long-range         = %-4s\n",dummy[0]);
  fprintf(fp,"tCNC__log_> lr-tol                 = %4.2f %4.2f\n",inp->lr_tol[0],inp->lr_tol[1]);
  fprintf(fp,"tCNC__log_> min-constr             = %d\n",inp->min);
  fprintf(fp,"tCNC__log_> max-constr             = %d\n",inp->max);
  fprintf(fp,"tCNC__log_> hbond-max-dist         = %4.2f\n",inp->hbond_max_dist);
  fprintf(fp,"tCNC__log_> hbond-min-angle        = %4.2f\n",inp->hbond_min_angle);
  for(i=0;i<3;i++){
    if(inp->use_hbond_protection[i]) {
      strcpy(dummy[i],"yes");
      dummy[i][3] = 0;
    }
    else {
      strcpy(dummy[i],"no");
      dummy[i][2] = 0;
    }
  }
  fprintf(fp,"tCNC__log_> use-hbond-protection   = %-4s %-4s %-4s\n",dummy[0],dummy[1],dummy[2]);
  fprintf(fp,"tCNC__log_> solvation-max          = %4.2f %4.2f %4.2f\n",
          inp->solvation_max[0],inp->solvation_max[1],inp->solvation_max[2]);
  fprintf(fp,"tCNC__log_> hphob-dist             = %4.2f\n",inp->hphob_dist);
  fprintf(fp,"\n");
}

/*==========================================================*/
void print_input_file(FILE *fp,t_input *inp)
{
  shstr dummy[3];
  int i;
  fprintf(fp,"; tCONCOORD PARAMETER FILE\n");
  fprintf(fp,";====================================\n");
  fprintf(fp,";\n");
  fprintf(fp,";\n");
  fprintf(fp,"; PARAMETERS\n");
  fprintf(fp,";=======================================\n");
  fprintf(fp,"\n");
  fprintf(fp,"; BOND FLEXIBILIY (A)\n");
  fprintf(fp,"\n");
  fprintf(fp,"bond-tol               = %5.2f\n",inp->bond_tol);
  fprintf(fp,"\n");
  fprintf(fp,"; ANGLE FLEXIBILIY (A)\n");
  fprintf(fp,"\n");
  fprintf(fp,"angle-tol              = %5.2f\n",inp->angle_tol);
  fprintf(fp,"\n");
  fprintf(fp,"; ANGLE FLEXIBILIY (deg)\n");
  fprintf(fp,"\n");
  fprintf(fp,"angle-sig              = %5.2f\n",inp->angle_sig);
  fprintf(fp,"\n");
  fprintf(fp,"; PLANARITY TOLERANCE\n");
  fprintf(fp,"\n");
  fprintf(fp,"plan-tol               = %6.4f\n",inp->plan_tol);
  fprintf(fp,"\n");
  fprintf(fp,"; BUMP TOLERANCE (A)\n");
  fprintf(fp,"\n");
  fprintf(fp,"bump-tol               = %5.2f\n",inp->bump_tol);
  fprintf(fp,"\n");
  fprintf(fp,"bump-14tol               = %5.2f\n",inp->bump_14tol);
  fprintf(fp,"\n");
  fprintf(fp,"; FLEXIBILITY FOR RESTRICTED DIHEDRALS (%%)\n");
  fprintf(fp,"\n");
  fprintf(fp,"restr-dihed-tol        = %5.2f\n",inp->restr_dihed_tol);
  fprintf(fp,"\n");
  fprintf(fp,"; FLEXIBILITY FOR NON-RESTRICTED DIHEDRALS\n");
  fprintf(fp,"\n");
  fprintf(fp,"non-restr-dihed-tol    = %5.2f\n",inp->non_restr_dihed_tol);
  fprintf(fp,"\n");
/*   fprintf(fp,"; USE PACKING EVALUATION\n");  */
/*   fprintf(fp,"; (keeps well-packed regions together)\n");  */
/*   fprintf(fp,"\n");  */
/*    if(inp->use_packing) strcpy(dummy[0],"yes");  */
/*    else strcpy(dummy[0],"no");  */
/*    fprintf(fp,"use-packing            = %5s\n",dummy[0]);  */
/*    fprintf(fp,"\n");  */
/*    fprintf(fp,"pack-limit             = %5.2f\n",inp->pack_limit);  */
/*   fprintf(fp,"; MINIMUM PACKING SCORE\n"); */
/*   fprintf(fp,"\n"); */
/*   fprintf(fp,"pack-limit             = %5.2f\n",inp->pack_limit); */
/*   fprintf(fp,"\n"); */
  if(inp->use_hbonds) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"use-hbonds                  = %5s\n",dummy[0]);
  if(inp->use_hbond_angles) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"use-hbond-angles            = %5s\n",dummy[0]);
  if(inp->use_sidechain) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"use-sidechain               = %5s\n",dummy[0]);
  if(inp->use_hydrophobics) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"use-hydrophobics            = %5s\n",dummy[0]);
  if(inp->use_network) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"use-network                 = %5s\n",dummy[0]);
  fprintf(fp,"network-tol                 = %8.3f %8.3f\n",inp->network_tol[0],inp->network_tol[1]);
  if(inp->use_packing) strcpy(dummy[0],"yes"); 
  else strcpy(dummy[0],"no"); 
  fprintf(fp,"use-packing                 = %5s\n",dummy[0]); 
  fprintf(fp,"\n"); 
  fprintf(fp,"pack-limit                  = %5.2f\n",inp->pack_limit); 

  if(inp->use_close_pairs) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"use-close-pairs             = %5s\n",dummy[0]);
  fprintf(fp,"close-pairs-fix-dist        = %8.3f\n",inp->close_pairs_fix_dist);
  if(inp->use_neighbor_ca) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"use-neighbor-ca             = %5s\n",dummy[0]);
  if(inp->use_long_range) strcpy(dummy[0],"yes");
  else strcpy(dummy[0],"no");
  fprintf(fp,"use-long-range              = %5s\n",dummy[0]);
  fprintf(fp,"lr-tol                      = %8.3f %8.3f\n",inp->lr_tol[0],inp->lr_tol[1]);
  fprintf(fp,"min-constr                  = %d\n",inp->min);
  fprintf(fp,"max-constr                  = %d\n",inp->max);


  fprintf(fp,"; MAX. DIST FOR HYDROGEN BONDS (A)\n");
  fprintf(fp,"\n");
  fprintf(fp,"hbond-max-dist         = %5.2f\n",inp->hbond_max_dist);
  fprintf(fp,"\n");
  fprintf(fp,"; MIN. ANGLE FOR HYDROGEN BONDS\n");
  fprintf(fp,"\n");
  fprintf(fp,"hbond-min-angle        = %5.2f\n",inp->hbond_min_angle);
  fprintf(fp,"\n");
  fprintf(fp,"; CHECK HYDROPHOBIC PROTECTION OF HYDROGEN\n");
  fprintf(fp,"; BONDS\n");
  fprintf(fp,"; backbone-backbone  sidechain-backbone   sidechain-sidechain\n");
  fprintf(fp,"\n");
  for(i=0;i<3;i++){
    if(inp->use_hbond_protection[i]) {
      strcpy(dummy[i],"yes");
      dummy[i][3] = 0;
    }
    else {
      strcpy(dummy[i],"no");
      dummy[i][2] = 0;
    }
  }
  fprintf(fp,"use-hbond-protection   = %5s %5s %5s\n",dummy[0],dummy[1],dummy[2]);
  fprintf(fp,"\n");
/*   fprintf(fp,"; CALCULATE SOLVATION SCORE USING ATOMS\n"); */
/*   fprintf(fp,"; WITHIN # ANGSTROEM AROUND THE HYDROGEN\n"); */
/*   fprintf(fp,"\n"); */
/*   fprintf(fp,"solvation-rad         = %5.2f %5.2f %5.2f\n", */
/*           inp->solvation_rad[0],inp->solvation_rad[1],inp->solvation_rad[2]); */
/*   fprintf(fp,"\n"); */
  fprintf(fp,"; MAXIMUM SOLVATION SCORE\n");
  fprintf(fp,"; (if the number is larger, the hydrogen bond\n");
  fprintf(fp,";  will be ignored )\n");
  fprintf(fp,"solvation-max         = %5.2f %5.2f %5.2f\n",
          inp->solvation_max[0],inp->solvation_max[1],inp->solvation_max[2]);
  fprintf(fp,"\n");
  fprintf(fp,"; HYDROPHOBIC CLUSTERS\n");
  fprintf(fp,"; (3 hydrophobic atoms from 3 different\n");
  fprintf(fp,";  residues form a cluster if\n");
  fprintf(fp,";  the distances a-b, a-c and b-c are smaller than #)\n");
  fprintf(fp,"\n");
  fprintf(fp,"hphob-dist             = %5.2f\n",inp->hphob_dist);
  fprintf(fp,"\n");
/*   fprintf(fp,"; RADIUS OF GYRATION\n"); */
/*   fprintf(fp,"Rg         = %5.2f %5.2f %5.2f\n", */
/*           inp->radius_of_gyration[XX], inp->radius_of_gyration[YY], inp->radius_of_gyration[ZZ]); */
/*   fprintf(fp,"Rg-tol         = %5.2f %5.2f %5.2f\n", */
/*           inp->radius_of_gyration_tol[XX], inp->radius_of_gyration_tol[YY], inp->radius_of_gyration_tol[ZZ]); */
  
  
  fclose(fp);
}
/*==========================================================*/
t_input *read_cnc_input(FILE *log,char *filename)
{
#define MAXPTR 254

  int ninp;
  int i,k;
  int nel;
  char *ptr[MAXPTR];
  char warn[STRLEN];
  char error[STRLEN];
  
  t_input *inp = input_init();
  if(fopen(filename,"r") == NULL){
    fprintf(stderr,"tCNC__log_> no parameter file given!! Using default values \n");
    sprintf(warn,"no paramater file given!! Using default values \n");
    CNCwarn(log,warn);
    sprintf(warn," -> THIS MIGHT NOT BE WHAT YOU WANT\n");      
    CNCwarn(log,warn);
    return inp;
  }
  t_inpfile *inf = read_inpfile(filename,&ninp);
  char dum[STRLEN][3];
  fprintf(stderr,"tCNC__log_> Reading input parameters: %s\n",filename);
  for(k=0;k<ninp;k++){
    if(strcmp(inf[k].name,"bond-tol") == 0)
      inp->bond_tol = atof(inf[k].value);
    else if(strcmp(inf[k].name,"angle-tol") == 0)
      inp->angle_tol = atof(inf[k].value);
    else if(strcmp(inf[k].name,"angle-sig") == 0)
      inp->angle_sig = atof(inf[k].value);
    else if(strcmp(inf[k].name,"plan-tol") == 0)
      inp->plan_tol = atof(inf[k].value);
    else if(strcmp(inf[k].name,"bump-tol") == 0)
      inp->bump_tol = atof(inf[k].value);
    else if(strcmp(inf[k].name,"bump-14tol") == 0)
      inp->bump_14tol = atof(inf[k].value);
    else if(strcmp(inf[k].name,"restr-dihed-tol") == 0)
      inp->restr_dihed_tol = atof(inf[k].value);
    else if(strcmp(inf[k].name,"non-restr-dihed-tol") == 0)
      inp->non_restr_dihed_tol = atof(inf[k].value);
    else if(strcmp(inf[k].name,"use-packing") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_packing = TRUE;
      else inp->use_packing = FALSE;
    else if(strcmp(inf[k].name,"pack-limit") == 0)
      inp->pack_limit = atof(inf[k].value);
    else if(strcmp(inf[k].name,"use-hbonds") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_hbonds = TRUE;
      else inp->use_hbonds = FALSE;
    else if(strcmp(inf[k].name,"use-hbond-angles") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_hbond_angles = TRUE;
      else inp->use_hbond_angles = FALSE;
    else if(strcmp(inf[k].name,"use-sidechain") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_sidechain = TRUE;
      else inp->use_hydrophobics = FALSE;
    else if(strcmp(inf[k].name,"use-hydrophobics") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_hydrophobics = TRUE;
      else inp->use_hydrophobics = FALSE;
    else if(strcmp(inf[k].name,"use-network") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_network = TRUE;
      else inp->use_network = FALSE;
    else if(strcmp(inf[k].name,"use-close-pairs") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_close_pairs = TRUE;
      else inp->use_close_pairs = FALSE;
    else if(strcmp(inf[k].name,"close-pairs-fix-dist") == 0)
      inp->close_pairs_fix_dist = atof(inf[k].value);
    else if(strcmp(inf[k].name,"use-neighbor-ca") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_neighbor_ca = TRUE;
      else inp->use_neighbor_ca = FALSE;
    else if(strcmp(inf[k].name,"use-long-range") == 0)
      if(strcmp(inf[k].value,"yes") == 0)
        inp->use_long_range = TRUE;
      else inp->use_long_range = FALSE;
    else if(strcmp(inf[k].name,"min-constr") == 0)
      inp->min = atoi(inf[k].value);
    else if(strcmp(inf[k].name,"max-constr") == 0)
      inp->max = atoi(inf[k].value);
    else if(strcmp(inf[k].name,"lr-tol") == 0) {
      nel = str_nelem(inf[k].value,MAXPTR,ptr);
      if(nel != 2){
        sprintf(error,"Invalid input: %s = %s\n",inf[k].name,inf[k].value);
        fatal_error(error);
      }
      for(i=0;i<nel;i++){
        inp->lr_tol[i] = atof(ptr[i]);
      }
    }

    else if(strcmp(inf[k].name,"network-tol") == 0) {
      nel = str_nelem(inf[k].value,MAXPTR,ptr);
      if(nel != 2){
        sprintf(error,"Invalid input: %s = %s\n",inf[k].name,inf[k].value);
        fatal_error(error);
      }
      for(i=0;i<nel;i++){
        inp->network_tol[i] = atof(ptr[i]);
      }
    }

    else if(strcmp(inf[k].name,"hbond-max-dist") == 0)
      inp->hbond_max_dist = atof(inf[k].value);
    else if(strcmp(inf[k].name,"hbond-min-angle") == 0)
      inp->hbond_min_angle = atof(inf[k].value);
    else if(strcmp(inf[k].name,"use-hbond-protection") == 0){
      nel = str_nelem(inf[k].value,MAXPTR,ptr);
      if(nel != 3){
        sprintf(error,"Invalid input: %s = %s\n",inf[k].name,inf[k].value);
        fatal_error(error);
      }
      for(i=0;i<nel;i++){
        if(strcmp(ptr[i],"yes") == 0) inp->use_hbond_protection[i] = TRUE;
        else  inp->use_hbond_protection[i] = FALSE;
      }
    }
    else if(strcmp(inf[k].name,"solvation-rad") == 0){
#ifdef GMX_DOUBLE
      sscanf(inf[k].value,"%lf %lf %lf",&inp->solvation_rad[0],
             &inp->solvation_rad[1],&inp->solvation_rad[2]);
#else
      sscanf(inf[k].value,"%f %f %f",&inp->solvation_rad[0],
             &inp->solvation_rad[1],&inp->solvation_rad[2]);
#endif
    }
#ifdef GMX_DOUBLE
    else if(strcmp(inf[k].name,"solvation-max") == 0){
      sscanf(inf[k].value,"%lf %lf %lf",&inp->solvation_max[0],
             &inp->solvation_max[1],&inp->solvation_max[2]);
    }
#else
    else if(strcmp(inf[k].name,"solvation-max") == 0){
      sscanf(inf[k].value,"%f %f %f",&inp->solvation_max[0],
             &inp->solvation_max[1],&inp->solvation_max[2]);
    }
#endif
    else if(strcmp(inf[k].name,"hphob-dist") == 0){
      inp->hphob_dist = atof(inf[k].value);
    }
    else{
      char estr[STRLEN];
      sprintf(estr,"Input Error: Unknown option %s\n",inf[k].name);
      CNCerr(log,estr);
    }
  }
  return inp;
}


/*==========================================================*/

