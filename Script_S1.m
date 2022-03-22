% Script to incorporate 1A01's GR-dependent biomass composition in FBA
% simulations (instead of using 1A01's biomass composition on glucose to
% build a canonical biomass reaction for all conditions). First, a growth
% rate on a given carbon source is guessed, and the corresponding biomass
% composition incorporated into the model via the biomass reaction. FBA
% then uses the model to predict an optimal growth rate. If this optimal
% growth rate matches the initially guessed growth rate, the script stops.
% Otherwise, the cycle repeats itself, with the FBA-predicted growth rate
% generating a new biomass composition (and, by extension, a new biomass
% reaction in the model), until the script converges (i.e. until the input
% and output growth rates match).

clear
close all
addpath('/Library/gurobi811/mac64/matlab') % Gurobi path
load model_1A01 % Load model

% Set up FBA
b = zeros(1,size(m.S,1)); % FBA parameter
contypes = char( '='*ones( 1,size(m.S,1) ) )'; % FBA parameter
constr_changes = []; % FBA parameter
bmIndex = find(strcmp('BiomassRxn',m.rxns)); % Identify biomass reaction in model
c = zeros(1,size(m.S,2)); % Define objective to optimize in FBA (namely, biomass production)
c(bmIndex) = 1;

% Identify carbon exchange reaction and set flux
cIndex = find(strcmp('EX_gal',m.rxns)); % Galactose exchange reaction
m.lb(cIndex) = -3.08; % Measured uptake flux on galactose is -3.08

% Guess growth rate
gr_in = 1;

%% BODY OF SCRIPT

gr_out = 0; % Initialize
while round(gr_in,2) ~= round(gr_out,2) % Repeat until input and output growth rates match
    % After first iteration, change output growth rate to input
    if gr_out ~= 0
        gr_in = gr_out;
    end
    
    % Display input growth rate
    disp(['Input g.r. is ' num2str(gr_in)])
    
    % Calculate protein and RNA biomass fractions based on OD (see
    % DM&R&P_1A01 file)
    protein_OD = -0.0928*gr_in + 0.3768; 
    RNA_OD = 0.0498*gr_in + 0.0547; 
    
    % Convert OD-fractions to DM-fractions
    DM = -0.135*gr_in + 0.63; 
    protein_DM = protein_OD/DM;
    RNA_DM = RNA_OD/DM;
    
    % Calculate biomass reaction based on protein and RNA DM-fractions
    biorxn = calc_biorxn(protein_DM,RNA_DM);
    
    % Edit model's biomass reaction
    m = edit_biorxn(biorxn,m);
    
    % Run FBA
    [~, gr_out] = FBA_Gurobi(c,m.S,b,contypes,m.lb,m.ub,constr_changes);
    
    % Display output growth rate
    disp(['Output g.r. is ' num2str(gr_out)]) 
end

%% FUNCTIONS

% Function to calculate biomass reaction based on GR-dependent protein/RNA fractions (see Supplementary Table 1 for full description)
function c = calc_biorxn(protein,RNA)

% Column C in Supplementary Table 1
osmolytes = 0.0586;
DNA = 0.0327*(1-RNA-protein-osmolytes)/0.2046;
murein = 0.0263*(1-RNA-protein-osmolytes)/0.2046;
LPS = 0.0358*(1-RNA-protein-osmolytes)/0.2046;
lipids = 0.0959*(1-RNA-protein-osmolytes)/0.2046;
inorganic_ions = 0.0105*(1-RNA-protein-osmolytes)/0.2046;
soluble_pools = 0.0034*(1-RNA-protein-osmolytes)/0.2046;
macro_perc = [protein DNA RNA murein LPS lipids inorganic_ions soluble_pools osmolytes];

% Column E in Supplementary Table 1
molmono_molmacro = [0.083 0.042 0.044 0.057 0.010 0.045 0.065 0.067 0.021 0.063 0.101 0.054 0.027 0.042 0.038 0.071 0.056 0.012 0.031 0.071 0.280 0.220 0.220 0.280 0.213 0.302 0.217 0.269 1 1 0.4590 0.5410 0.7143 0.0476 0.0317 0.0190 0.0286 0.0286 0.0127 0.0127 0.0127 0.0127 0.0127 0.0190 0.0159 0.0159 0.000576 0.001831 0.000447 0.000223 0.000223 0.000223 0.000223 0.000223 0.000223 0.000223 0.000223 0.000223 0.000223 0.000055 0.000223 0.000223 0.000223 0.102389078 0.897610922];

% Columns G to L in Supplementary Table 1
C = [3	6	4	4	3	5	5	2	6	6	6	6	5	9	5	3	4	11	9	5	10	9	10	10	9	10	9	10	77	84	37	37	0	0	0	0	0	0	0	0	0	0	0	0	0	0	21	21	21	27	19	20	20	20	12	49	8	34	30	55	10	15	17	5	5	0	0];
H = [7	15	8	6	7	10	8	5	9	13	13	15	11	11	9	7	9	12	11	11	12	10	12	13	12	12	11	12	117	148	74	70	0	4	0	0	0	0	0	0	0	0	0	0	0	1	32	26	25	31	21	21	24	21	16	76	8	30	27	89	8	23	20	10	8	2	1];
N = [1	4	2	1	1	2	1	1	3	1	1	2	1	1	1	1	1	2	1	1	5	3	5	2	3	5	2	5	15	2	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	7	7	7	9	7	7	7	7	4	0	1	4	3	0	0	6	4	2	1	0	0];
O = [2	2	3	4	2	3	4	2	2	2	2	2	2	2	2	3	3	2	3	2	12	13	13	14	14	14	15	13	40	37	8	8	0	0	0	0	0	0	0	0	4	0	0	0	4	4	16	14	17	15	6	6	6	7	7	4	6	4	15	7	6	5	6	3	4	1	7];
P = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3	3	3	3	3	3	3	3	0	2	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	3	2	3	2	0	0	0	0	2	0	1	0	0	2	0	0	0	0	0	0	2];
S = [0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	0	0];

% Column M in Supplementary Table 1
for i = 1 : length(C)
    if i == 33
        MW(i) = 38.9637;
    elseif i == 35
        MW(i) = 23.985;
    elseif i == 36
        MW(i) = 39.9626;
    elseif i == 37
        MW(i) = 55.9349;
    elseif i == 38
        MW(i) = 55.9349;
    elseif i == 39
        MW(i) = 63.546;
    elseif i == 40
        MW(i) = 54.938;
    elseif i == 41
        MW(i) = 159.94;
    elseif i == 42
        MW(i) = 58.9332;
    elseif i == 43
        MW(i) = 63.9291;
    elseif i == 44
        MW(i) = 34.9689;
    else
        MW(i) = (C(i)*12.011)+(H(i)*1.008)+(O(i)*15.999)+(N(i)*14.007)+(P(i)*30.974)+(S(i)*32.066);
    end
end

% Column O in Supplementary Table 1
for i = 1 : length(MW)-2
    if i < 21
        MW_corr(i) = MW(i) - MW(end-1);
    elseif i > 20 && i < 29
        MW_corr(i) = MW(i) - MW(end);
    else
        MW_corr(i) = MW(i);
    end
end

% Column P in Supplementary Table 1
for i = 1 : length(MW_corr)
    if i < 47 || i == 64 || i == 65
        gmono_molmacro(i) = molmono_molmacro(i)*MW_corr(i);
    end
end

% Column Q in Supplementary Table 1
gmacro_molmacro(1) = sum(gmono_molmacro(1:20));
gmacro_molmacro(2) = sum(gmono_molmacro(21:24));
gmacro_molmacro(3) = sum(gmono_molmacro(25:28));
gmacro_molmacro(4) = gmono_molmacro(29);
gmacro_molmacro(5) = gmono_molmacro(30);
gmacro_molmacro(6) = sum(gmono_molmacro(31:32));
gmacro_molmacro(7) = sum(gmono_molmacro(33:46));
gmacro_molmacro(8) = sum(gmono_molmacro(64:65));

% Column R in Supplementary Table 1
for i = 1 : length(gmono_molmacro)
    if i < 21
        gmono_gmacro(i) = gmono_molmacro(i)/gmacro_molmacro(1);
    elseif i > 20 && i < 25
        gmono_gmacro(i) = gmono_molmacro(i)/gmacro_molmacro(2);
    elseif i > 24 && i < 29
        gmono_gmacro(i) = gmono_molmacro(i)/gmacro_molmacro(3);
    elseif i == 29
        gmono_gmacro(i) = gmono_molmacro(i)/gmacro_molmacro(4);
    elseif i == 30
        gmono_gmacro(i) = gmono_molmacro(i)/gmacro_molmacro(5);
    elseif i > 30 && i < 33
        gmono_gmacro(i) = gmono_molmacro(i)/gmacro_molmacro(6);
    elseif i > 32 && i < 47
        gmono_gmacro(i) = gmono_molmacro(i)/gmacro_molmacro(7);
    elseif i > 63 && i < 66
        gmono_gmacro(i) = gmono_molmacro(i)/gmacro_molmacro(8);
    end
end

% Column S in Supplementary Table 1
for i = 1 : length(gmono_gmacro)
    if i < 21
        gmono_gCDW(i) = gmono_gmacro(i)*macro_perc(1);
    elseif i > 20 && i < 25
        gmono_gCDW(i) = gmono_gmacro(i)*macro_perc(2);
    elseif i > 24 && i < 29
        gmono_gCDW(i) = gmono_gmacro(i)*macro_perc(3);
    elseif i == 29
        gmono_gCDW(i) = gmono_gmacro(i)*macro_perc(4);
    elseif i == 30
        gmono_gCDW(i) = gmono_gmacro(i)*macro_perc(5);
    elseif i > 30 && i < 33
        gmono_gCDW(i) = gmono_gmacro(i)*macro_perc(6);
    elseif i > 32 && i < 47
        gmono_gCDW(i) = gmono_gmacro(i)*macro_perc(7);
    elseif i > 63 && i < 66
        gmono_gCDW(i) = gmono_gmacro(i)*macro_perc(9);
    end
end

% Column T in Supplementary Table 1
for i = 1 : length(gmono_gmacro)
    if i < 47
        mmolmono_gCDW(i) = gmono_gCDW(i)/MW_corr(i)*1000;
    elseif i > 46 && i < 64
        mmolmono_gCDW(i) = molmono_molmacro(i)*macro_perc(8)/0.0034;
    elseif i > 63 && i < 66
        mmolmono_gCDW(i) = gmono_gCDW(i)/MW_corr(i)*1000;
    end
end
mmolmono_gCDW(70) = 15.81;
mmolmono_gCDW(71) = 15.81;
mmolmono_gCDW(72) = sum(mmolmono_gCDW(21:24));
mmolmono_gCDW(73) = sum(mmolmono_gCDW(25:28));
mmolmono_gCDW(74) = sum(mmolmono_gCDW(1:20));
mmolmono_gCDW(66) = mmolmono_gCDW(70) + mmolmono_gCDW(28);
mmolmono_gCDW(67) = mmolmono_gCDW(64) + mmolmono_gCDW(6);
mmolmono_gCDW(68) = mmolmono_gCDW(65) + mmolmono_gCDW(7);
mmolmono_gCDW(69) = mmolmono_gCDW(70) - mmolmono_gCDW(74);
mmolmono_gCDW(75) = mmolmono_gCDW(73) + mmolmono_gCDW(72);
mmolmono_gCDW(76) = mmolmono_gCDW(70) - mmolmono_gCDW(46);

% Column U in Supplementary Table 1
c(1:5) = -mmolmono_gCDW(1:5);
c(6:25) = -mmolmono_gCDW(8:27);
c(26:27) = -mmolmono_gCDW(29:30);
c(28) = -(0.0223/(0.0223+0.0415))*mmolmono_gCDW(31);
c(29) = -(0.0415/(0.0223+0.0415))*mmolmono_gCDW(31);
c(30) = -(0.0263/(0.0263+0.0489))*mmolmono_gCDW(32);
c(31) = -(0.0489/(0.0263+0.0489))*mmolmono_gCDW(32);
c(32:44) = -mmolmono_gCDW(33:45);
c(45:61) = -mmolmono_gCDW(47:63);
c(62:65) = -mmolmono_gCDW(66:69);
c(66:67) = mmolmono_gCDW(70:71);
c(68:69) = mmolmono_gCDW(75:76);
c = c';

end

% Function to edit model's biomass reaction
function m_out = edit_biorxn(c,m)

% Column V in Supplementary Table 1
bio_mets = {'L-ALPHA-ALANINE[c]' 'ARG[c]' 'ASN[c]' 'L-ASPARTATE[c]' 'CYS[c]'...
    'GLY[c]' 'HIS[c]' 'ILE[c]' 'LEU[c]' 'LYS[c]' 'MET[c]' 'PHE[c]' 'PRO[c]'...
    'SER[c]' 'THR[c]' 'TRP[c]' 'TYR[c]' 'VAL[c]' 'DATP[c]' 'DCTP[c]' 'DGTP[c]'...
    'TTP[c]' 'CTP[c]' 'GTP[c]' 'UTP[c]' 'CPD0-2278[p]' 'KDO2-LIPID-IVA[p]' 'CPD-12819[c]'...
    'CPD-12819[p]' 'CPD-17086[c]' 'CPD-17086[p]' 'K+[c]' 'AMMONIUM[c]' 'MG+2[c]'...
    'CA+2[c]' 'FE+2[c]' 'FE+3[c]' 'CU+2[c]' 'MN+2[c]' 'CPD-3[c]' 'CO+2[c]' 'ZN+2[c]'...
    'CL-[c]' 'SULFATE[c]' 'CO-A[c]' 'NAD[c]' 'NADP[c]' 'FAD[c]' 'THF-GLU-N[c]'...
    'METHYLENE-THF-GLU-N[c]' '5-METHYL-THF-GLU-N[c]' '10-FORMYL-DIHYDROFOLATE-GLU-N[c]'...
    'THIAMINE-PYROPHOSPHATE[c]' 'CPD-9956[c]' 'PYRIDOXAL_PHOSPHATE[c]' 'PROTOHEME[c]'...
    'ENTEROBACTIN[c]' 'UNDECAPRENYL-DIPHOSPHATE[c]' 'CHORISMATE[c]' 'S-ADENOSYLMETHIONINE[c]'...
    'RIBOFLAVIN[c]' 'ATP[c]' 'GLN[c]' 'GLT[c]' 'WATER[c]' 'ADP[c]' 'PROTON[c]' 'PPI[c]' 'Pi[c]'};
bmIndex = find(strcmp('BiomassRxn',m.rxns)); % Identify biomass reaction in model
m.S(:,bmIndex) = 0; % Delete model's existing biomass reaction
for i = 1 : length(bio_mets)
    metIndex = find(strcmp(bio_mets{i},m.mets));
    if isempty(metIndex)
        error('Metabolite not found')
    end
    m.S(metIndex, bmIndex) = c(i); % Assign new biomass reaction coefficients
end
m_out = m;

end
