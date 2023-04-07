%% 
% To reproduce figures from Ziegler et al.
% 1. Please replace basedir with'<DATA_DOWNLOAD_LOCATION>\Ziegler_et_al_Supplementary_Data\Data_tables'
% 2. Please add the helper_functions directory to the Matlab path.
% 3. Please download and add the violinplot libary to the Matlab path: https://github.com/bastibe/Violinplot-Matlab

clear 
close all

%Ziegler_et_al_in_vivo_ephys\Data_tables
basedir='\\lsdf02.urz.uni-heidelberg.de\sd19b001\PainData\Analysis_PainPaper1\Ziegler_et_al_Supplementary_Data\Data_tables'
basedir='C:\Users\Rebecca Mease\Documents\local_sds_working\Ziegler_et_al_Supplementary_Data\Data_tables'
cd(basedir);
%% Figure 1
clearvars -except basedir
PlotFig1B;
%Figure 1c in supplementary data file: Fig1c_S1_light_response_overlay_m380_S1_1205.fig

clearvars -except basedir
%Figure 1d (and related Figure S8a, Figure S8b)
VPL_S1_variable_powers_script

%% Figure 2
% All data in Supplemental Excel file 

%% Figure 3
clearvars -except basedir
%Figure 3b-h (and related Figure S8c1-c3,d)
VPL_ML_script;
%Fix violin for RP

%% Figure 4
clearvars -except basedir
%Figure 4b-f (and related Figure S8e)
S1_ML_L6_script
%% Figure 5
clearvars -except basedir
%Figure 5c,d 
S1_stgt_L5_L_script
%Remaining data in Supplemental Excel file 

%% Figure 6
%Figure 6b
PlotFig6B;
%Figure 6c in supplementary data file: Fig6c_L5_l5_chr_light_response_m393pt1_snipped_514.fig

%% Figure 7
% All data in Supplemental Excel file 
%% Figure S1
% All data in figure.
%% Figure S2
%Figure S2 a,b
PlotFigS2;
%Figure S2 c
PlotFigS2c;
%% Figure S3-6
% All data in Supplemental Excel file 

%% Figure S7
% All data in figure.

%% Figure S8
%Figure S8a, Figure S8b: see section above for Figure 1.
%Figure S8c1-c3, Figure S8d: see section above for Figure 3.
%Figure S8e, see section above for Figure 4.

%% Figure S9
% All data in Supplemental Excel file 
%% Figure S10
% All data in Supplemental Excel file 
%% Figure S11
% All data in Supplemental Excel file 
%% Figure S12 POm data
%Figure S12b, Figure S12c, Figure S12d
clearvars -except basedir
POm_ML_script
