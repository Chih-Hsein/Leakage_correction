function C = extended_tofts_model(Ktrans, ve, Cp, Tvec)
% Compute leaked tissue concentration defined in Extended Tofts Model.
%
%   Inputs:
%       Ktrans : Transfer constant from plasma to EES
%       ve     : Extravascular, extracellular volume fraction
%       Cp     : Plasma concentration (AIF) over time [vector]
%       Tvec   : Corresponding time vector [vector]
%
%   Output:
%       C : Estimated tissue concentration as a function of time.
%
%   This function computes:
%
%   C(t) = Ktrans * ∫ [Cp(τ) * exp(- Ktrans/ve * (t - τ)) dτ ]
%
%   Implementation uses trapezoidal integration (trapz).
%

    % Ensure column vectors
    Cp   = Cp(:);
    Tvec = Tvec(:);

    % Preallocate for speed
    C = zeros(1, numel(Tvec));

    for k = 1:numel(Tvec)
        % Current time up to index k
        T_k  = Tvec(1:k);
        Cp_k = Cp(1:k);

        % Exponential weighting
        F = Cp_k .* exp((-Ktrans ./ ve) .* (T_k(end) - T_k));

        % Perform numerical integration using trapz
        if numel(T_k) == 1
            % If there's only one point, trapz won't integrate properly
            % because it interprets single values differently
            M = 0;
        else
            M = trapz(T_k, F);
        end

        C(k) = Ktrans .* M;
    end

    C = C(:);

end
