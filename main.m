%--------------------------------------------------------------------------
% DEMO: DCE-to-DSC Leakage Correction in Brain Tumor Imaging
%--------------------------------------------------------------------------
% This MATLAB script demonstrates a region-based signal ratio model fitting 
% approach for both DCE and DSC data. It then uses estimated parameters 
% to compensate for contrast leakage in DSC images.
%
% Example data are taken from a patient with glioma.
% © 2024 Chih-Hsien Tseng <c.tseng@tudelft.nl>
%
%--------------------------------------------------------------------------

clc;
clear all
close all;

%% USER-DEFINED PARAMETERS
% These would normally be set according to your actual scanning parameters.

% DCE imaging parameters
T1_initial    = 1.98;    % [s] measured tissue T1
FA_DCE        = 25;      % [degrees] flip angle for DCE
TR_DCE        = 0.0027;  % [s] repetition time for DCE
temp_res      = 2;       % [s] temporal resolution for both DCE and DSC

% DSC imaging parameters
TR_DSC        = 2;       % [s] repetition time for DSC
TE_DSC        = 0.045;   % [s] echo time for DSC
FA_DSC        = 90;      % [degrees] flip angle for DSC

% Relaxivity for gadolinium
r1            = 4.5;     % [mM^-1 s^-1] longitudinal relaxivity
% The code also uses DSCfitP.r2t for transverse relaxivity in tissue.

%% SECTION 1: Load Measured Data (Example)

% Measured DCE signal ratio curve in enhancing tumor
SigR_DCE = [
    1.08,0.99,0.96,0.99,1.01,1.00,0.99,0.99,1.02,0.98,...
    0.98,1.00,1.01,1.19,1.60,1.95,1.95,1.90,1.74,1.66,...
    1.65,1.68,1.70,1.67,1.70,1.70,1.73,1.71,1.74,1.75,...
    1.76,1.78,1.80,1.80,1.80,1.80,1.79,1.85,1.85,1.84,...
    1.83,1.86,1.86,1.83,1.85,1.84,1.84,1.89,1.88,1.89,...
    1.91,1.89,1.95,1.92,1.93,1.91,1.91,1.94,1.97,1.92,...
    1.94,1.96,1.94,1.96,1.96,1.97,1.95,2.00,2.00,2.02,...
    1.98,2.03,1.98,2.02,1.99,2.00,2.07,2.01,2.01,2.01,...
    2.02,2.01,2.07,2.04,2.02,2.05,2.06,2.05,2.07,2.12,...
    2.05,2.10,2.13,2.08,2.10,2.10,2.14,2.13,2.12,2.19
];

% Measured AIF from DCE, after inflow and partial volume correction
% AIF correction following the method in https://doi.org/10.1002/nbm.5225
AIF = [
    0, -0.04, -0.02, -0.01, 0.03, 0.03, 0.02, 0.04, -0.01, -0.02,...
    -0.02, 0, 0.03, 2.14, 6.22, 8.56, 6.90, 4.54, 2.89, 1.87,...
    1.59, 1.68, 1.88, 2.04, 2.05, 1.94, 1.90, 1.73, 1.76, 1.76,...
    1.83, 1.79, 1.75, 1.79, 1.69, 1.61, 1.55, 1.48, 1.57, 1.55,...
    1.47, 1.46, 1.40, 1.42, 1.38, 1.33, 1.32, 1.29, 1.28, 1.31,...
    1.35, 1.39, 1.36, 1.35, 1.26, 1.19, 1.17, 1.17, 1.21, 1.18,...
    1.17, 1.13, 1.12, 1.16, 1.14, 1.20, 1.15, 1.16, 1.13, 1.08,...
    1.13, 1.14, 1.18, 1.12, 1.09, 1.13, 1.09, 1.11, 1.03, 1.05,...
    1.03, 1.03, 1.05, 1.01, 0.96, 1.03, 1.00, 1.01, 1.04, 1.03,...
    1.08, 1.02, 1.04, 1.05, 0.95, 0.92, 1.00, 1.01, 1.06, 1.01
];

% Measured DSC signal ratio curve in enhancing tumor (manually time-aligned)
SigR_DSC = [
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.01, 1.00, 1.01, ...
    1.00, 1.00, 0.98, 0.88, 0.72, 0.62, 0.61, 0.65, 0.70, 0.76, ...
    0.81, 0.84, 0.85, 0.85, 0.86, 0.85, 0.85, 0.86, 0.85, 0.85, ...
    0.85, 0.86, 0.86, 0.87, 0.86, 0.87, 0.87, 0.87, 0.86, 0.86, ...
    0.86, 0.86, 0.86, 0.87
];

%% SECTION 2: Define Time Vectors
% Ensure both are truncated to match their respective signal lengths.

numDCE = length(SigR_DCE);
t      = 0 : temp_res : temp_res*(numDCE-1);

numDSC = length(SigR_DSC);
t_DSC  = t(1:numDSC);  % Reuse the first portion of t to match DSC length

AIF_DSC = AIF(1:numDSC);  % Align with DSC

%% SECTION 3: DCE Signal Ratio Model Fit
disp('=== Fitting DCE signal ratio model ===');

% Fit the defined model to the measured DCE signal ratio
[DCEfitP, ~, ~, fitC_DCE] = fit_sigR_DCE(SigR_DCE, AIF, t, TR_DCE, FA_DCE, T1_initial);

% Compute the leaked tissue concentration from the Extended Tofts Model
LeakC = extended_tofts_model(DCEfitP.ktrans, DCEfitP.ve, AIF, t);

%% SECTION 4: DSC Signal Ratio Model Fit
disp('=== Fitting DSC signal ratio model with r2* estimation ===');

% Fit the DSC ratio model, re-using the DCE-derived ktrans, ve, vc
[DSCfitP, ~, ~, fitC_DSC] = fit_sigR_DSC( ...
    SigR_DSC, AIF_DSC, t_DSC, TR_DSC, TE_DSC, FA_DSC, ...
    DCEfitP.ktrans, DCEfitP.ve, DCEfitP.vc);

%% SECTION 5: Compute and Plot Leakage Correction
disp('=== Computing and plotting the leakage correction ===');

% Truncate leaked concentration to match DSC length
LeakC_DSC = LeakC(1:numDSC)';

% Convert DSC signal ratio to Delta R2*
DeltaR2st = (-1 / TE_DSC) .* log(SigR_DSC);

% The T1 leakage term
T1_leakTerm = (1 / TE_DSC) .* log( ...
    (1 - exp(-TR_DSC / DSCfitP.T10t) .* exp(-TR_DSC * r1 * LeakC_DSC)) ./ ...
    (1 - exp(-TR_DSC / DSCfitP.T10t)) );

% The T2* leakage term
T2star_leakTerm = DSCfitP.r2t .* LeakC_DSC;

% Corrected ΔR2*: eq. (DeltaR2* + T1term - T2*term)
DeltaR2st_corr = DeltaR2st + T1_leakTerm - T2star_leakTerm;

% --- Plot the contaminated ΔR2* and leakage components ---
figure('Name','Leakage Terms');
hold on; grid on;
plot(t_DSC, DeltaR2st,       'LineWidth',1.5, 'DisplayName','Measured \DeltaR_2^*');
plot(t_DSC, T1_leakTerm,     'LineWidth',1.5, 'DisplayName','T_1 Leakage Term');
plot(t_DSC, T2star_leakTerm, 'LineWidth',1.5, 'DisplayName','T_2^* Leakage Term');
xlabel('Time (s)');
ylabel('\Delta R_2^* (1/s)');
legend('Location','Best');
title('Measured \DeltaR_2^* vs. Estimated Leakage Components');

% --- Plot the final corrected ΔR2* curve ---
figure('Name','Corrected \DeltaR_2^*');
hold on; grid on;
plot(t_DSC, DeltaR2st,       'b-o','LineWidth',1.5, 'DisplayName','Measured \DeltaR_2^*');
plot(t_DSC, DeltaR2st_corr,  'r-o','LineWidth',1.5, 'DisplayName','Corrected \DeltaR_2^*');
xlabel('Time (s)');
ylabel('\Delta R_2^* (1/s)');
legend('Location','Best');
title('Leakage Correction in DSC');

disp('=== Processing complete. ===');