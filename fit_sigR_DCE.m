function [FittingParameters, gof, output, fitcurve] = fit_sigR_DCE(SigR, Cp, t, Tr, FA, T1)
%Fit a DCE signal ratio curve based on the proposed signal ratio model with seperated blood and tissue component.
%
%   Inputs:
%       SigR : Measured DCE signal ratio (vector)
%       Cp   : Arterial/blood concentration curve (vector)
%       t    : Time vector corresponding to Cp (same length as Cp)
%       Tr   : Repetition time for DCE sequence [s]
%       FA   : Flip angle for DCE sequence [degrees]
%       T1   : Measured T1 of the tissue [s]
%
%   Outputs:
%       FittingParameters : Fit result object with estimated parameters 
%                           (ktrans, vc, ve).
%       gof              : Goodness-of-fit metrics (see 'fit' doc).
%       output           : Structure with additional fitting info.
%       fitcurve         : Model-predicted DCE signal ratio using the fitted parameters.
%
%   FittingParameters :
%       Ktrans : Transfer constant from plasma to EES
%       vc     : Blood plasma volume fraction
%       ve     : Extravascular, extracellular volume fraction
%
%   See also: extended_tofts_model

    % Ensure column vectors
    Cp = Cp(:);
    t  = t(:);

    % T1 relaxivity for gadolinium-based agent 
    r1 = 4.5;           

    % Inverse T1 times for blood & tissue
    T10b_re = 1 / 1.8;  % Blood T1 is assumed ~1.8 s => reciprocal
    T10t_re = 1 / T1;   % Reciprocal of measured tissue T1

    % Tissue concentration according to ETM
    c = @(Cp, ktrans, ve, t) ...
        extended_tofts_model(ktrans, ve, Cp, t);

    % Exponential factor for T1-based signal
    E1 = @(T1_re) exp(-Tr .* T1_re);

    % Blood T1 relaxation
    % (0.55 = 1 - Hct, with Hct=0.45)
    T1_blood = @(Cp) (T10b_re + r1 .* 0.55 .* Cp);

    % Tissue T1 relaxation 
    T1_tissue = @(Cp, ktrans, ve, t) ...
        (T10t_re + r1 .* c(Cp, ktrans, ve, t));

    % Signal model: weighted sum of blood + tissue
    Sig_model = @(vc, T1b_re, T1t_re) ...
        ( vc      .* (1 - E1(T1b_re)) ./ (1 - cosd(FA).*E1(T1b_re)) + ...
         (1 - vc) .* (1 - E1(T1t_re)) ./ (1 - cosd(FA).*E1(T1t_re)) );
    
    % Signal ratio model
    D = @(ktrans, vc, ve, Cp, t) ...
        Sig_model( ...
            vc, ...
            T1_blood(Cp), ...
            T1_tissue(Cp, ktrans, ve, t) ) ...
        ./ ...
        Sig_model(vc, T10b_re, T10t_re);

    % Define bounds and initial guesses
    LB = [0,     0,  0];
    UB = [0.01,  1,  1];
    st = [0.001, 0.01, 0.01];  % [ktrans, vc, ve]

    % Fit options
    fitOpts = fitoptions('Method','NonlinearLeastSquares',...
                         'Lower',LB,'Upper',UB,'StartPoint',st);

    % Create fit type:
    f = fittype(D, ...
        'independent',  {'Cp','t'}, ...
        'coefficients', {'ktrans','vc','ve'}, ...
        'options',       fitOpts);

    % Perform the fit
    [FittingParameters, gof, output] = fit([Cp, t], SigR', f);

    % Display the fitted parameters
    disp(FittingParameters);

    % Generate the fitted curve 
    fitcurve = D(FittingParameters.ktrans, FittingParameters.vc, ...
                 FittingParameters.ve, Cp, t);

end
