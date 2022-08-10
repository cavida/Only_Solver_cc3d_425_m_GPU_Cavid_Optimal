import antimony
import os

def antimonyToSBML(ant):
    """ Convert Antimony to SBML string.

    :param ant: Antimony string or file
    :type ant: str | file
    :return: SBML
    :rtype: str
    """
    antimony.clearPreviousLoads()
    antimony.freeAll()
    try:
        isfile = os.path.isfile(ant)
    except ValueError:
        isfile = False
    if isfile:
        code = antimony.loadAntimonyFile(ant)
    else:
        code = antimony.loadAntimonyString(ant)
#     _checkAntimonyReturnCode(code)
    if code < 0:
        raise Exception('Cavid: After converting Antimony to SBML, code <0! Conversion problem')
    mid = antimony.getMainModuleName()
    return antimony.getSBMLString(mid)


simple_model='''
A -> B; 0.000000001*A
B -> C; 1000000*B 
A=1
'''

Mukti_model_antimony='''
// Created by libAntimony v2.12.0.3
model *EndMT_v1()

  // Compartments and Species:
  compartment Medium, Cytoplasm, Nucleus;
  species $SD093 in Cytoplasm, TGFBeta in Medium, TGFBeta_Nonspecific in Medium;
  species TGFBeta_Cytoplasm in Cytoplasm, TGFBR_Complex in Cytoplasm, TGFBR_Complex_Cytoplasm in Cytoplasm;
  species TGFBR1_Surface in Cytoplasm, TGFBR1_Cytoplasm in Cytoplasm, TGFBR1_Surface_Inhibited in Cytoplasm;
  species TGFBR2_Surface in Cytoplasm, TGFBR2_Cytoplasm in Cytoplasm, Smad2 in Cytoplasm;
  species Smad2_Nucleus in Nucleus, pSmad2 in Cytoplasm, pSmad2_Nucleus in Nucleus;
  species Smad4 in Cytoplasm, Smad4_Nucleus in Nucleus, pSmad24 in Cytoplasm;
  species pSmad24_Nucleus in Nucleus, pSmad22 in Cytoplasm, pSmad22_Nucleus in Nucleus;
  species Smad3 in Cytoplasm, Smad3_Nucleus in Nucleus, pSmad3 in Cytoplasm;
  species pSmad3_Nucleus in Nucleus, pSmad34 in Cytoplasm, pSmad34_Nucleus in Nucleus;
  species pSmad33 in Cytoplasm, pSmad33_Nucleus in Nucleus, $Snail in Nucleus;
  species $alphaSMA in Nucleus, $Slug in Nucleus, $PECAM1 in Nucleus, $MMP9 in Nucleus, $init_PECAM1 in Nucleus, last_alphaSMA in Nucleus; // Cavid added init_PECAM1 and last_alphaSMA
  species $mRNASnail in Nucleus, $mRNAalphaSMA in Nucleus, $mRNASlug in Nucleus;
  species $mRNAPECAM1 in Nucleus, $mRNAMMP9 in Nucleus, $Amino_Acids_Placeholder in Cytoplasm;
  species $Degraded_Proteins in Cytoplasm;

  // Assignment Rules:
  r_xj_Snail := SNAI1_gene*v_x_max*(mRNASnail/(mRNASnail + k_sat_protein));
  r_xj_alphaSMA := ACTA2_gene*v_x_max*(mRNAalphaSMA/(mRNAalphaSMA + k_sat_protein));
  r_xj_Slug := Slug_gene*v_x_max*(mRNASlug/(mRNASlug + k_sat_protein));
  r_xj_PECAM1 := PECAM1_gene*v_x_max*(mRNAPECAM1/(mRNAPECAM1 + k_sat_protein));
  r_xj_MMP9 := MMP9_gene*v_x_max*(mRNAMMP9/(mRNAMMP9 + k_sat_protein));
  u_pSmad24_Snail := f_pSmad24_Snail;
  r_mj := v_t_max*(number_of_genes/(number_of_genes + k_sat_mrna));
  u_Snail_alphaSMA := f_Snail_alphaSMA;
  u_pSmad34_Slug := f_pSmad34_Slug;
  u_Slug_PECAM1 := 1 - f_Slug_PECAM1;
  u_Slug_MMP9 := f_Slug_MMP9;
  totalNuclearPSmads := pSmad2_Nucleus + pSmad24_Nucleus + pSmad3_Nucleus + pSmad34_Nucleus + 2*pSmad22_Nucleus + 2*pSmad33_Nucleus;
  f_pSmad24_Snail := k_control_activation_pSmad24^hill_coeff_pSmad24*pSmad24_Nucleus^hill_coeff_pSmad24/(1 + k_control_activation_pSmad24^hill_coeff_pSmad24*pSmad24_Nucleus^hill_coeff_pSmad24);
  f_Snail_alphaSMA := k_control_activation_Snail^hill_coeff_Snail*Snail^hill_coeff_Snail/(1 + k_control_activation_Snail^hill_coeff_Snail*Snail^hill_coeff_Snail);
  f_pSmad34_Slug := k_control_activation_pSmad34^hill_coeff_pSmad34*pSmad34_Nucleus^hill_coeff_pSmad34/(1 + k_control_activation_pSmad34^hill_coeff_pSmad34*pSmad34_Nucleus^hill_coeff_pSmad34);
  f_Slug_MMP9 := k_control_activation_Slug_MMP9^hill_coeff_Slug_MMP9*Slug^hill_coeff_Slug_MMP9/(1 + k_control_activation_Slug_MMP9^hill_coeff_Slug_MMP9*Slug^hill_coeff_Slug_MMP9);
  f_Slug_PECAM1 := k_control_activation_Slug_PECAM1^hill_coeff_Slug_PECAM1*Slug^hill_coeff_Slug_PECAM1/(1 + k_control_activation_Slug_PECAM1^hill_coeff_Slug_PECAM1*Slug^hill_coeff_Slug_PECAM1);

  // Rate Rules:
  Snail' = times(r_xj_Snail - Snail*(mu_growth_rate + theta_xj));
  alphaSMA' = times(r_xj_alphaSMA - alphaSMA*(mu_growth_rate + theta_xj));
  Slug' = times(r_xj_Slug - Slug*(mu_growth_rate + theta_xj));
  PECAM1' = times(r_xj_PECAM1 - PECAM1*(mu_growth_rate + theta_xj));
  MMP9' = times(r_xj_MMP9 - MMP9*(mu_growth_rate + theta_xj));
  mRNASnail' = times(u_pSmad24_Snail*(SNAI1_gene*r_mj - mRNASnail*(mu_growth_rate + theta_mj))/1);
  mRNAalphaSMA' = times(u_Snail_alphaSMA*(ACTA2_gene*r_mj - mRNAalphaSMA*(mu_growth_rate + theta_mj))/1);
  mRNASlug' = times(u_pSmad34_Slug*(Slug_gene*r_mj - mRNASlug*(mu_growth_rate + theta_mj))/1);
  mRNAPECAM1' = 1*(u_Slug_PECAM1*(PECAM1_gene*r_mj - mRNAPECAM1*(mu_growth_rate + theta_mj))/1);
  mRNAMMP9' = times(u_Slug_MMP9*(MMP9_gene*r_mj - mRNAMMP9*(mu_growth_rate + theta_mj))/1);

  // Reactions:
  R1: TGFBeta -> TGFBeta_Nonspecific; Medium*(k_on_nonspecific*TGFBeta - k_off_nonspecific*TGFBeta_Nonspecific);
  R2: TGFBeta_Cytoplasm => $Degraded_Proteins; Cytoplasm*k_deg_TGFBeta*TGFBeta_Cytoplasm;
  R3: $Amino_Acids_Placeholder => TGFBR1_Surface; Cytoplasm*k_prod_TGFBR1;
  R4: $Amino_Acids_Placeholder => TGFBR2_Surface; Cytoplasm*k_prod_TGFBR2;
  R5: TGFBR1_Cytoplasm => TGFBR1_Surface; Cytoplasm*k_r*TGFBR1_Cytoplasm;
  R6: TGFBR2_Cytoplasm => TGFBR2_Surface; Cytoplasm*k_r*TGFBR2_Cytoplasm;
  R7: TGFBR1_Surface => TGFBR1_Cytoplasm; Cytoplasm*k_i*TGFBR1_Surface;
  R8: TGFBR2_Surface => TGFBR2_Cytoplasm; Cytoplasm*k_i*TGFBR2_Surface;
  R9: TGFBR1_Cytoplasm => $Degraded_Proteins; Cytoplasm*k_deg_TGFBR1*TGFBR1_Cytoplasm;
  R10: TGFBR2_Cytoplasm => $Degraded_Proteins; Cytoplasm*k_deg_TGFBR2*TGFBR2_Cytoplasm;
  R11: TGFBR1_Surface + TGFBR2_Surface + TGFBeta => TGFBR_Complex; Cytoplasm*ka_TGFBR_Complex*TGFBR1_Surface*TGFBR2_Surface*TGFBeta;
  R12: TGFBR_Complex => TGFBR_Complex_Cytoplasm; Cytoplasm*k_i*TGFBR_Complex;
  R13: TGFBR_Complex_Cytoplasm => $Degraded_Proteins; Cytoplasm*k_deg_TGFBR_Complex*TGFBR_Complex_Cytoplasm;
  R14: TGFBR_Complex_Cytoplasm => TGFBR1_Cytoplasm + TGFBR2_Cytoplasm + TGFBeta_Cytoplasm; Cytoplasm*k_diss_TGFBR_Complex*TGFBR_Complex_Cytoplasm;
  R15: TGFBR_Complex => $Degraded_Proteins; Cytoplasm*k_deg_TGFBR_Complex_Negative_Feedback*TGFBR_Complex*totalNuclearPSmads;
  R16: TGFBR_Complex_Cytoplasm + Smad2 => pSmad2 + TGFBR_Complex_Cytoplasm; Cytoplasm*k_phos*Smad2*TGFBR_Complex_Cytoplasm;
  R17: Smad4 + pSmad2 -> pSmad24; Cytoplasm*k_on*pSmad2*Smad4 - Cytoplasm*k_off*pSmad24;
  R18: Smad2 => Smad2_Nucleus; Cytoplasm*k_import_Smad2*Smad2;
  R19: Smad2_Nucleus => Smad2; Nucleus*k_export_Smad2*Smad2_Nucleus;
  R20: Smad4 => Smad4_Nucleus; Cytoplasm*k_import_Smad4*Smad4;
  R21: Smad4_Nucleus => Smad4; Nucleus*k_export_Smad4*Smad4_Nucleus;
  R22: pSmad2 => pSmad2_Nucleus; Cytoplasm*k_import_Smad2*pSmad2;
  R23: pSmad2_Nucleus => pSmad2; Nucleus*k_export_Smad2*pSmad2_Nucleus;
  R24: pSmad24 => pSmad24_Nucleus; Cytoplasm*k_import_Smad24*pSmad24;
  R25: pSmad24_Nucleus -> pSmad2_Nucleus + Smad4_Nucleus; Nucleus*(k_off*pSmad24_Nucleus - k_on*pSmad2_Nucleus*Smad4_Nucleus);
  R26: pSmad2_Nucleus => Smad2_Nucleus; Nucleus*k_dephos*pSmad2_Nucleus;
  R27: 2 pSmad2 -> pSmad22; Cytoplasm*(k_on*pSmad2^2 - k_off*pSmad22);
  R28: pSmad22 => pSmad22_Nucleus; Cytoplasm*k_import_Smad24*pSmad22;
  R29: pSmad22_Nucleus -> 2 pSmad2_Nucleus; Nucleus*(k_off*pSmad22_Nucleus - k_on*pSmad2_Nucleus^2);
  R30: TGFBR_Complex_Cytoplasm + Smad3 => pSmad3 + TGFBR_Complex_Cytoplasm; Cytoplasm*k_phos*Smad3*TGFBR_Complex_Cytoplasm;
  R31: Smad4 + pSmad3 -> pSmad34; Cytoplasm*k_on*pSmad3*Smad4 - Cytoplasm*k_off*pSmad34;
  R32: Smad3 => Smad3_Nucleus; Cytoplasm*k_import_Smad3*Smad3;
  R33: Smad3_Nucleus => Smad3; Nucleus*k_export_Smad3*Smad3_Nucleus;
  R34: pSmad3 => pSmad3_Nucleus; Cytoplasm*k_import_Smad3*pSmad3;
  R35: pSmad3_Nucleus => pSmad3; Nucleus*k_export_Smad3*pSmad3_Nucleus;
  R36: pSmad34 => pSmad34_Nucleus; Cytoplasm*k_import_Smad34*pSmad34;
  R37: pSmad34_Nucleus -> pSmad3_Nucleus + Smad4_Nucleus; Nucleus*(k_off*pSmad34_Nucleus - k_on*pSmad3_Nucleus*Smad4_Nucleus);
  R38: pSmad3_Nucleus => Smad3_Nucleus; Nucleus*k_dephos*pSmad3_Nucleus;
  R39: 2 pSmad3 -> pSmad33; Cytoplasm*(k_on*pSmad3^2 - k_off*pSmad33);
  R40: pSmad33 => pSmad33_Nucleus; Cytoplasm*k_import_Smad34*pSmad33;
  R41: pSmad33_Nucleus -> 2 pSmad3_Nucleus; Nucleus*(k_off*pSmad33_Nucleus - k_on*pSmad3_Nucleus^2);
  R42: TGFBR1_Surface + $SD093 -> TGFBR1_Surface_Inhibited; Cytoplasm*(k_on_sd093*TGFBR1_Surface*SD093 - k_off_sd093*TGFBR1_Surface_Inhibited);


// this section is commented out
//    u_pSmad24_Snail = (k_control_activation_pSmad24^hill_coeff_pSmad24)*(pSmad24_Nucleus^hill_coeff_pSmad24);
//    u_pSmad34_Slug = (k_control_activation_pSmad34^hill_coeff_pSmad34)*(pSmad34_Nucleus^hill_coeff_pSmad34);
//    u_Snail_alphaSMA = (k_control_activation_Snail^hill_coeff_Snail)*(Snail^hill_coeff_Snail);
//    u_Slug_MMP9 = (k_control_activation_Slug_MMP9^hill_coeff_Slug_MMP9)*(Slug^hill_coeff_Slug_MMP9);
//    u_Slug_PECAM1 = (k_control_activation_Slug_PECAM1^hill_coeff_Slug_PECAM1)*(Slug^hill_coeff_Slug_PECAM1);


  // Events:
  SD_093_Input: at time == 3600, t0=false: SD093 = 1000;
  SNAI1_gene_activation: at pSmad24_Nucleus >= tau_SNAI1, t0=false: SNAI1_gene = 1;
  SNAI1_gene_deactivation: at pSmad24_Nucleus < tau_SNAI1, t0=false: SNAI1_gene = 0;
  ACTA2_gene_activation: at Snail > tau_ACTA2, t0=false: ACTA2_gene = 1;
  ACTA2_gene_deactivation: at Snail <= tau_ACTA2, t0=false: ACTA2_gene = 0;
  Slug_gene_activation: at pSmad34_Nucleus > tau_Slug, t0=false: Slug_gene = 1;
  Slug_gene_deactivation: at pSmad34_Nucleus <= tau_Slug, t0=false: Slug_gene = 0;
  MMP9_gene_activation: at Slug > tau_MMP9, t0=false: MMP9_gene = 1;
  MMP9_gene_deactivation: at Slug <= tau_MMP9, t0=false: MMP9_gene = 0;
  PECAM1_gene_activation: at Slug > tau_PECAM1, t0=false: PECAM1_gene = 0;
  PECAM1_gene_deactivation: at Slug <= tau_PECAM1, t0=false: PECAM1_gene = 1;

  // Species initializations:
  SD093 = 0;
  TGFBeta = 0.0042 //TGFBeta = 0.0022;                                 
  TGFBeta_Nonspecific = 0;
  TGFBeta_Cytoplasm = 0;
  TGFBR_Complex = 0;
  TGFBR_Complex_Cytoplasm = 0;
  TGFBR1_Surface = 0.702;
  TGFBR1_Cytoplasm = 6.523;
  TGFBR1_Surface_Inhibited = 0;
  TGFBR2_Surface = 0.201;
  TGFBR2_Cytoplasm = 1.439;
  Smad2 = 30.3;
  Smad2_Nucleus = 14.25;
  pSmad2 = 0;
  pSmad2_Nucleus = 0;
  Smad4 = 50.8;
  Smad4_Nucleus = 50.8;
  pSmad24 = 0;
  pSmad24_Nucleus = 0;
  pSmad22 = 0;
  pSmad22_Nucleus = 0;
  Smad3 = 30.3;
  Smad3_Nucleus = 14.25;
  pSmad3 = 0;
  pSmad3_Nucleus = 0;
  pSmad34 = 0;
  pSmad34_Nucleus = 0;
  pSmad33 = 0;
  pSmad33_Nucleus = 0;
  Snail=0;
  alphaSMA = 0;
  Slug = 0;
  PECAM1 = 16.6;
  MMP9 = 0;
  mRNASnail = 0;
  mRNAalphaSMA = 0;
  mRNASlug = 0;
  mRNAPECAM1 = 4.981;
  mRNAMMP9 = 0;
  Amino_Acids_Placeholder = 0;
  Degraded_Proteins = 0;
  
// Cavid Additions; Also defined them in the compartment part
init_PECAM1 = 16.6;
last_alphaSMA = 0.0;
  

  // Compartment initializations:
  Medium = 7e-09;  // 2e-09 Previous value for the medium in Mukti's model
  Medium has volume;
  Cytoplasm = 2.27e-12;
  Cytoplasm has volume;
  Nucleus = 1e-12;
  Nucleus has volume;

  // Variable initializations:
  mu_growth_rate = 4.1805e-06;
  mu_growth_rate has ps;
  theta_xj = 0.000501;
  theta_xj has ps;
  SNAI1_gene = 0; 
  theta_mj = 4.1805e-05;
  theta_mj has ps;
  ACTA2_gene = 0;
  Slug_gene = 0;
  PECAM1_gene = 1;
  MMP9_gene = 0;
  tau_SNAI1 = 9.96008;
  tau_SNAI1 has nM;
  tau_ACTA2 = 0.00374;
  tau_ACTA2 has nM;
  tau_Slug = 9.96008;
  tau_Slug has nM;
  tau_MMP9 = 0.001;
  tau_MMP9 has nM;
  tau_PECAM1 = 0.001;
  tau_PECAM1 has nM;
  k_prod_TGFBR1 = 0.0002783;
  k_prod_TGFBR1 has nMp;
  k_deg_TGFBR1 = 4.267e-05;
  k_deg_TGFBR1 has ps;
  k_prod_TGFBR2 = 0.000316;
  k_prod_TGFBR2 has nMp;
  k_deg_TGFBR2 = 0.00022;
  k_deg_TGFBR2 has ps;
  k_r = 0.000555;
  k_r has ps;
  k_i = 0.00555;
  k_i has ps;
  ka_TGFBR_Complex = 1.9649;
  ka_TGFBR_Complex has third_order;
  k_deg_TGFBR_Complex = 4.267e-05;
  k_deg_TGFBR_Complex has ps;
  k_diss_TGFBR_Complex = 0.0007301;
  k_diss_TGFBR_Complex has ps;
  k_deg_TGFBR_Complex_Negative_Feedback = 0.0003894;
  k_deg_TGFBR_Complex_Negative_Feedback has ps;
  k_on_nonspecific = 0.000842355;
  k_on_nonspecific has ps;
  k_off_nonspecific = 0.033884;
  k_off_nonspecific has ps;
  k_deg_TGFBeta = 0.005783;
  k_deg_TGFBeta has ps;
  k_on_sd093 = 5;
  k_on_sd093 has pnMp;
  k_off_sd093 = 100;
  k_off_sd093 has ps;
  k_phos = 0.0008;
  k_phos has pnMp;
  k_dephos = 0.00656;
  k_dephos has ps;
  k_off = 0.016;
  k_off has ps;
  k_on = 0.003307;
  k_on has pnMp;
  k_import_Smad2 = 0.0026;
  k_import_Smad2 has ps;
  k_export_Smad2 = 0.01271;
  k_export_Smad2 has ps;
  k_import_Smad3 = 0.0026;
  k_import_Smad3 has ps;
  k_export_Smad3 = 0.01271;
  k_export_Smad3 has ps;
  k_import_Smad4 = 0.0026;
  k_import_Smad4 has ps;
  k_export_Smad4 = 0.005983;
  k_export_Smad4 has ps;
  k_import_Smad24 = 0.01481;
  k_import_Smad24 has ps;
  k_import_Smad34 = 0.01481;
  k_import_Smad34 has ps;
  totalNuclearPSmads has nM;
  k_sat_mrna = 3.37;
  k_sat_protein = 73.2;
  v_t_max = 0.00651;
  v_x_max = 0.01;
  number_of_genes = 2;
  k_control_activation_pSmad24 = 0.05;
  hill_coeff_pSmad24 = 3;
  k_control_activation_Snail = 0.05;
  hill_coeff_Snail = 1;
  k_control_activation_pSmad34 = 0.05;
  hill_coeff_pSmad34 = 3;
  k_control_activation_Slug_MMP9 = 0.05;
  hill_coeff_Slug_MMP9 = 0.1;
  k_control_activation_Slug_PECAM1 = 1.5;
  hill_coeff_Slug_PECAM1 = 1;

  // Other declarations:
  var r_xj_Snail, r_xj_alphaSMA, r_xj_Slug, r_xj_PECAM1, r_xj_MMP9, u_pSmad24_Snail;
  var SNAI1_gene, r_mj, u_Snail_alphaSMA, ACTA2_gene, u_pSmad34_Slug, Slug_gene;
  var u_Slug_PECAM1, PECAM1_gene, u_Slug_MMP9, MMP9_gene, totalNuclearPSmads;
  var f_pSmad24_Snail, k_control_activation_pSmad24, f_Snail_alphaSMA, f_pSmad34_Slug;
  var f_Slug_MMP9, f_Slug_PECAM1;
  const Medium, Cytoplasm, Nucleus, mu_growth_rate, theta_xj, theta_mj, tau_SNAI1;
  const tau_ACTA2, tau_Slug, tau_MMP9, tau_PECAM1, k_prod_TGFBR1, k_deg_TGFBR1;
  const k_prod_TGFBR2, k_deg_TGFBR2, k_r, k_i, ka_TGFBR_Complex, k_deg_TGFBR_Complex;
  const k_diss_TGFBR_Complex, k_deg_TGFBR_Complex_Negative_Feedback, k_on_nonspecific;
  const k_off_nonspecific, k_deg_TGFBeta, k_on_sd093, k_off_sd093, k_phos;
  const k_dephos, k_off, k_on, k_import_Smad2, k_export_Smad2, k_import_Smad3;
  const k_export_Smad3, k_import_Smad4, k_export_Smad4, k_import_Smad24, k_import_Smad34;
  const k_sat_mrna, k_sat_protein, v_t_max, v_x_max, number_of_genes, hill_coeff_pSmad24;
  const k_control_activation_Snail, hill_coeff_Snail, k_control_activation_pSmad34;
  const hill_coeff_pSmad34, k_control_activation_Slug_MMP9, hill_coeff_Slug_MMP9;
  const k_control_activation_Slug_PECAM1, hill_coeff_Slug_PECAM1;

  // Unit definitions:
  unit substance = 1e-9 mole;
  unit volume = litre;
  unit time_unit = time_unit;
  unit area = metre^2;
  unit length = metre;
  unit nM = 1e-9 mole / litre;
  unit ps = 1 / second;
  unit pnMp = litre / (1e-9 mole * second);
  unit nMp = 1e-9 mole / (second * litre);
  unit molecule = 1e-24 mole;
  unit third_order = 1 / ((1e-9 mole)^2 * second);
  unit extent = substance;

  // Display Names:
  time_unit is "time";
  ps is "persecond";
  pnMp is "pernMpersecond";
  nMp is "nMpersecond";
  third_order is "third_order_rate_constant";
  r_xj_Snail is "specific_rate_of_translation_for_Snail";
  mu_growth_rate is "specific_growth_rate";
  theta_xj is "degradation_constant_protein";
  r_xj_alphaSMA is "specific_rate_of_translation_for_alphaSMA";
  r_xj_Slug is "specific_rate_of_translation_for_Slug";
  r_xj_PECAM1 is "specific_rate_of_translation_for_PECAM1";
  r_xj_MMP9 is "specific_rate_of_translation_for_MMP9";
  u_pSmad24_Snail is "control_term_for_Snail_transcription";
  r_mj is "specific_rate_of_transcription";
  theta_mj is "degradation_constant_mRNA";
  u_Snail_alphaSMA is "control_term_for_ACTA2_transcription";
  u_pSmad34_Slug is "control_term_for_Slug_transcription";
  u_Slug_PECAM1 is "control_term_for_PECAM1_transcription";
  u_Slug_MMP9 is "control_term_for_MMP9_transcription";
  tau_SNAI1 is "Snail_activation_by_pSmad24_threshold";
  tau_ACTA2 is "ACTA2_activation_by_Snail_threshold";
  tau_Slug is "Slug_activation_by_pSmad34_threshold";
  tau_MMP9 is "MMP9_activation_by_Slug_threshold";
  tau_PECAM1 is "PECAM1_activation_by_Slug_threshold";
  k_prod_TGFBR1 is "k_prod_rate_of_TGFBR1_receptor_production";
  k_deg_TGFBR1 is "k_deg_rate_of_TGFBR1_receptor_degradation";
  k_prod_TGFBR2 is "k_prod_rate_of_TGFBR2_receptor_production";
  k_deg_TGFBR2 is "k_deg_rate_of_TGFBR2_receptor_degradation";
  k_r is "k_recycling_rate_for_receptors";
  k_i is "k_internalization_rate_for_receptors";
  k_phos is "k_phosphorylation_rate";
  k_dephos is "k_dephosphorylation_rate";
  k_off is "k_off_Smad_complex_dissociation_rate";
  k_on is "k_on_Smad_complex_formation_rate";
  k_import_Smad2 is "k_Smad2_nuclear_import";
  k_export_Smad2 is "k_Smad2_nuclear_export";
  k_import_Smad3 is "k_Smad3_nuclear_import";
  k_export_Smad3 is "k_Smad3_nuclear_export";
  k_import_Smad4 is "k_Smad4_nuclear_import";
  k_export_Smad4 is "k_Smad4_nuclear_export";
  k_import_Smad24 is "k_Smad24_complex_import";
  k_import_Smad34 is "k_Smad34_complex_import";
  totalNuclearPSmads is "total_Smads_in_Nucleus";
  k_sat_mrna is "saturation_constant_mRNA";
  k_sat_protein is "saturation_constant_protein";
  v_t_max is "maximum_gene_expression_rate";
  v_x_max is "maximum_translation_rate";
  f_pSmad24_Snail is "transfer_function_pSmad24_Snail";
  k_control_activation_pSmad24 is "species_gain_parameter_pSmad24";
  hill_coeff_pSmad24 is "cooperativity_pSmad24";
  f_Snail_alphaSMA is "transfer_function_Snail_ACTA2";
  k_control_activation_Snail is "species_gain_parameter_Snail";
  hill_coeff_Snail is "cooperativity_Snail";
  f_pSmad34_Slug is "transfer_function_pSmad34_Slug";
  k_control_activation_pSmad34 is "species_gain_parameter_pSmad34";
  hill_coeff_pSmad34 is "cooperativity_pSmad34";
  f_Slug_MMP9 is "transfer_function_Slug_MMP9";
  k_control_activation_Slug_MMP9 is "species_gain_parameter_Slug_MMP9";
  hill_coeff_Slug_MMP9 is "cooperativity_Slug_MMP9";
  f_Slug_PECAM1 is "transfer_function_Slug_PECAM1";
  k_control_activation_Slug_PECAM1 is "species_gain_parameter_Slug_PECAM1";
  hill_coeff_Slug_PECAM1 is "cooperativity_Slug_PECAM1";
  R1 is "Reaction_1_TGFB_NonSpecific_Binding";
  R2 is "Reaction_2_TGFB_Degradation_by_Endocytosis";
  R3 is "Reaction_3_TGFBR1_Production";
  R4 is "Reaction_4_TGFBR2_Production";
  R5 is "Reaction_5_TGFBR1_Recycling";
  R6 is "Reaction_6_TGFBR2_Recycling";
  R7 is "Reaction_7_TGFBR1_Internalization_through_Endocytosis";
  R8 is "Reaction_8_TGFBR2_Internalization_through_Endocytosis";
  R9 is "Reaction_9_TGFBR1_Degradation";
  R10 is "Reaction_10_TGFBR2_Degradation";
  R11 is "Reaction_11_TGFBR_Complex_Activation";
  R12 is "Reaction_12_TGFBR_Complex_Internalization";
  R13 is "Reaction_13_TGFBR_Complex_Degradation";
  R14 is "Reaction_14_TGFBR_Complex_Dissociation";
  R15 is "Reaction_15_TGFBR_Complex_Degradation_due_to_Negative_Feedback_of_Smads";
  R16 is "Reaction_16_Smad2_phosphorylation";
  R17 is "Reaction_17_Smad2_Smad4_complex_formation";
  R18 is "Reaction_18_Smad2_Nuclear_Import";
  R19 is "Reaction_19_Smad2_Nuclear_Export";
  R20 is "Reaction_20_Smad4_Nuclear_Import";
  R21 is "Reaction_21_Smad4_Nuclear_Export";
  R22 is "Reaction_22_pSmad2_Nuclear_Import";
  R23 is "Reaction_23_pSmad2_Nuclear_Export";
  R24 is "Reaction_24_pSmad24_Nuclear_Import";
  R25 is "Reaction_25_pSmad24_Dissociation";
  R26 is "Reaction_26_pSmad2_Dephosphorylation";
  R27 is "Reaction_27_pSmad2_Dimer_Formation";
  R28 is "Reaction_28_pSmad22_Nuclear_Import";
  R29 is "Reaction_29_pSmad2_Dimer_Dissociation";
  R30 is "Reaction_30_Smad3_phosphorylation";
  R31 is "Reaction_31_Smad3_Smad4_complex_formation";
  R32 is "Reaction_32_Smad3_Nuclear_Import";
  R33 is "Reaction_33_Smad3_Nuclear_Export";
  R34 is "Reaction_34_pSmad3_Nuclear_Import";
  R35 is "Reaction_35_pSmad3_Nuclear_Export";
  R36 is "Reaction_36_pSmad34_Nuclear_Import";
  R37 is "Reaction_37_pSmad34_Dissociation";
  R38 is "Reaction_38_pSmad3_Dephosphorylation";
  R39 is "Reaction_39_pSmad3_Dimer_Formation";
  R40 is "Reaction_40_pSmad33_Nuclear_Import";
  R41 is "Reaction_41_pSmad3_Dimer_Dissociation";
  R42 is "Reaction_42_TGFBR1_Inhibition";
end

EndMT_v1 is "TGFB_Smad_Signaling"

'''
