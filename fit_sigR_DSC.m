function [fitp, gof, output, fitresult] = fit_sigR_DSC( ...
    SigR, Cp, t, Tr, Te, FA, ktrans, ve, vc)
% Fit DSC signal ratio curve, solving for T10t and r2t.
%
%   [fitp, gof, output, fitresult] = fit_sigR_DSC( ...
%       SigR, Cp, t, Tr, Te, FA, ktrans, ve, vc )
%
%   Inputs:
%       SigR   : DSC signal ratio (vector, length N)
%       Cp     : Arterial/blood concentration curve (N)
%       t      : Time vector (N)
%       Tr     : Repetition time [s]
%       Te     : Echo time [s]
%       FA     : Flip angle [deg]
%       ktrans : from DCE fit
%       ve     : from DCE fit
%       vc     : from DCE fit
%
%   This function fits T10t, r2t as new coefficients. The others (ktrans, ve, vc)
%   are treated as known constants from the prior DCE step.
%
%   Outputs:
%       fitp     : struct with T10t, r2t
%       gof      : goodness-of-fit metrics
%       output   : additional info
%       fitresult: fitted DSC signal ratio curve
%
%   See also: extended_tofts_model

    % Convert to column vectors
    SigR = SigR(:);
    Cp   = Cp(:);
    t    = t(:);

    % Example constants
    refC      = mean(Cp(end-3:end)); % AIF tail average
    r1        = 4.5; % T1 relaxivity for blood and tissue 
    r2Blood   = 6; % T2* relaxivity for blood 
    T10b_re   = 1 / 1.8 + 4.5*0.55*refC; % Taking remaining concentration in blood into account 
    T20b_re   = 1 / 0.02; % Assumed initial T2* value in blood: 20ms
    T20t_re   = 1 / 0.02; % Assumed initial T2* value in tissue: 20ms

    % Tissue concentration 
    c = @ (Cp, t) extended_tofts_model(ktrans, ve, Cp, t);

    E1 = @(T1_re) exp(-Tr .* T1_re);
    E2 = @(T2_re) exp(-Te .* T2_re);

    % Blood T1, T2*
    T1_blood = @(Cp) T10b_re + r1 .* 0.55 .* Cp;
    T2_blood = @(Cp) T20b_re + r2Blood .* 0.55 .* Cp;

    % Tissue T1, T2* with unknown T10t, r2t
    T1_tissue = @(T10t, Cp, t) ...
                (1./T10t) + r1.*c(Cp, t);
    T2_tissue = @(r2t, Cp, t) ...
                T20t_re + r2t.*(c(Cp, t) + vc.*0.55.*Cp);

    % DSC signal model
    Sig_model = @(T1b_re, T2b_re, T1t_re, T2t_re) ...
        ( vc .* (1 - E1(T1b_re)) ./ (1 - cosd(FA).*E1(T1b_re)) .* E2(T2b_re) + ...
         (1 - vc).* (1 - E1(T1t_re)) ./ (1 - cosd(FA).*E1(T1t_re)) .* E2(T2t_re) );

    % The key: fit T10t, r2t; use Cp, t as independent
    D = @(T10t, r2t, Cp, t) ...
        Sig_model( ...
            T1_blood(Cp), ...
            T2_blood(Cp), ...
            T1_tissue(T10t, Cp, t), ...
            T2_tissue(r2t,  Cp, t) ) ...
        ./ ...
        Sig_model(T10b_re, T20b_re, 1./T10t, T20t_re);

    % Fit bounds for T10t and r2t
    LB = [0,   0];
    UB = [5, 300];
    st = [0.5, 30];

    fitOpts = fitoptions('Method','NonlinearLeastSquares',...
                         'Lower',LB, 'Upper',UB, 'StartPoint',st);

    f = fittype(D, ...
        'coefficients', {'T10t','r2t'}, ...
        'independent',  {'Cp','t'}, ...
        'options',      fitOpts);

    % Perform the fit
    [fitp, gof, output] = fit([Cp, t], SigR, f);

    % The fitted curve
    fitresult = D(fitp.T10t, fitp.r2t, Cp, t);

    disp(fitp);

end
