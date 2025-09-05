function [q, q05, q95, Q, Q05, Q95, s_m10_HF] = q_deridderEq24(Hm0, Hm0_HF, Tm10_HF, Rc, varargin)
% Q_DERIDDER_EQ24_HF  Mean overtopping discharge (De Ridder et al., 2024, Eq. 24)
%                     using HF parameters to compute s_{m-1,0,HF} internally,
%                     with log-normal 5–95% predictive bands. This formula was calibrated with a short–long (SS–IG) band cutoff of 0.45/T_{m−1,0}
%
% INPUTS
%   Hm0       [m] : significant wave height at the toe (full spectrum)
%   Hm0_HF    [m] : significant wave height of the HF part of the spectrum
%   Tm10_HF   [s] : mean spectral period T_{m-1,0} of the HF part
%   Rc        [m] : freeboard
%
% NAME-VALUE (optional)
%   'gamma_f'   [-] : roughness factor (default = 0.55 for rock, calibrated by De Ridder)
%   'sigma_ln'  [-] : ln-scatter (RMSLE proxy) for 5–95% bands (default = 0.59)
%
% OUTPUTS
%   q      [m^3/s/m] : mean overtopping discharge (dimensional)
%   q05,q95           : 5% / 95% predictive bands (dimensional)
%   Q      [-]        : non-dimensional discharge (q / sqrt(g*Hm0^3))
%   Q05,Q95           : 5% / 95% bands (non-dimensional)
%   s_m10_HF [-]      : computed HF steepness = (2π/g)*Hm0_HF / Tm10_HF^2
%
% NOTES
%   - s_{m-1,0,HF} = (2π/g) * Hm0_HF / Tm10_HF^2  (equiv. Hm0_HF / L0_HF with L0_HF = g T^2 / 2π)
%   - Bands use a log-normal envelope: Q05=Q*exp(-z*σ_ln), Q95=Q*exp(+z*σ_ln), z=1.645 (5–95%).
%
% Example:
%   [q,q05,q95,Q,Q05,Q95,sHF] = q_deridderEq24_HF(3.0, 1.2, 10.0, 2.0, ...
%       'gamma_f',0.55, 'sigma_ln',0.64);
%
% Author: César Esparza A., 5/9/2025 University of Edinburgh
% -------------------------------------------------------------------------

    % Defaults
    g         = 9.81;
    gamma_f   = 0.55;     % set to 0.55 for rubble-mound if desired
    sigma_ln  = 0.59;    % optimistic RMSLE from De Ridder

    % Parse name-value pairs
    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'gamma_f'
                gamma_f = varargin{k+1};
            case 'sigma_ln'
                sigma_ln = varargin{k+1};
            case 'placement'
                placement = lower(varargin{k+1});
        end
    end

    % Compute HF steepness from HF parameters
    s_m10_HF = (2*pi/g) .* (Hm0_HF ./ (Tm10_HF.^2));   % s_{m-1,0,HF}

    % Core Eq. (24) - choose placement of s^{0.32}
            Q = 0.74 .* exp( -8.51 .* ( Rc ./ (Hm0 .* gamma_f) ) .* (s_m10_HF.^0.32) );


    % Log-normal 5–95% predictive bands
    z = 1.645;
    fac_low  = exp(-z * sigma_ln);
    fac_high = exp( z * sigma_ln);
    Q05 = fac_low  .* Q;
    Q95 = fac_high .* Q;

    % Convert to dimensional discharge using full-spectrum Hm0
    q   = Q   .* sqrt(g .* Hm0.^3);
    q05 = Q05 .* sqrt(g .* Hm0.^3);
    q95 = Q95 .* sqrt(g .* Hm0.^3);
end
