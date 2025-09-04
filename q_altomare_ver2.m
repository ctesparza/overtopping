function [q, q_5, q_95, Q, Q_5, Q_95] = q_altomare(Htoe, Tm, Rc, htoe, alpha, m)
% Q_ALTOMARE  Mean wave overtopping discharge on a sloping structure using
%             Altomare et al. (2016) "equivalent slope" concept, with
%             one–sigma scatter expressed as 5% and 95% confidence curves.
%
% ───────────────────────────────────────────────────────────────────────────
% PURPOSE 
%   Implements the Altomare et al. (2016) formulation for the *dimensionless*
%   mean overtopping discharge, Q = q / sqrt(g * Htoe^3), on a smooth
%   impermeable slope under *normal* wave attack, accounting for a foreshore
%   (m) via an *equivalent slope* (tan δ) computed iteratively from run-up.
%
%   The code assumes all reduction factors γ (roughness, obliquity, porosity,
%   breaking) are equal to 1. If your structure is a rubble mound, waves are
%   oblique, or other effects apply, adjust γ-factors following EurOtop (2018).
%
% INPUTS (SI units unless noted)
%   Htoe  [m]   : Spectral significant wave height at the structure toe.
%   Tm    [s]   : Mean spectral period at the toe, T_{m-1,0}. (Ensure this
%                 truly is T_{m-1,0}; if estimated, e.g. via Hofland et al. 2017,
%                 do that pre-processing upstream.)
%   Rc    [m]   : Freeboard (vertical distance from SWL to crest).
%   htoe  [m]   : Water depth at the toe (relative to SWL; positive if submerged).
%   alpha [-]   : Structure slope in V:H units (e.g., 1/2 for 1V:2H).
%   m     [-]   : Foreshore slope in V:H units (e.g., 1/30 for 1V:30H).
%
% OUTPUTS
%   q     [m^3/s/m] : Mean overtopping discharge (dimensional).
%   q_5   [m^3/s/m] : Lower confidence (5% curve; conservative low).
%   q_95  [m^3/s/m] : Upper confidence (95% curve; conservative high).
%   Q     [-]       : Dimensionless mean discharge, q / sqrt(g * Htoe^3).
%   Q_5   [-]       : Dimensionless 5% curve.
%   Q_95  [-]       : Dimensionless 95% curve.
%
% METHOD (brief)
%   1) Compute deep-water wavelength L0 ≈ 1.56 * Tm^2 (linear theory).
%   2) Iterate run-up Ru using the equivalent slope tan δ:
%        - Build a composite horizontal length Lslope over foreshore + structure.
%        - Update tan δ = (1.5 Htoe + Ru) / Lslope.
%        - Compute Iribarren with equivalent slope, ξ_δ = tan δ / sqrt(Htoe/L0).
%        - Update Ru = (4.0 - 1.5 / sqrt(ξ_δ)) * γ_break * Htoe.
%      Iterate until relative change in Ru falls below 1% (with |·|).
%   3) With converged ξ_δ, evaluate
%        Q = 10^c * exp( - Rc / ( Htoe * (0.33 + 0.022 * ξ_δ) ) ),
%      with c = -0.791 (Altomare 2016 regression).
%   4) Confidence envelopes are modelled by shifting the intercept:
%        c_± = c ± z * σ_c, with σ_c = 0.294 and z = 1.64 (≈95% one-sided).
%      Then compute Q_5, Q_95 and convert to dimensional q_5, q_95.
%
% NOTES & CAVEATS
%   - Valid for smooth slopes and normal attack unless γ-factors are applied.
%   - The term (1.5*Htoe - htoe) in Lslope should remain physically meaningful;
%     ensure your geometry is consistent with Altomare et al. (2016) assumptions.
%   - L0 = 1.56 Tm^2 is a deep-water approximation; here it only affects ξ_δ.
%
% REFERENCES
%   - Altomare, C. et al. (2016) " Wave overtopping of sea dikes with very shallow foreshores". 
%   - EurOtop (2018) "Manual on wave overtopping of sea defences and related 
%     structures". Reduction factors γ_f, γ_β, γ_nu, γ_break.
%
% Author: César Esparza A. University of Edinburgh (rewritten, commented, and minor fixes)
%

%% Constants & γ-factors (set to 1 here; adjust if needed)
g           = 9.81;   % [m/s^2]
gamma_f     = 1.0;    %  % roughness (not used here explicitly)
gamma_beta  = 1.0;    % % obliquity (not used here explicitly)
gamma_nu    = 1.0;    %  % permeability/porosity (not used here)
gamma_break = 1.0;               % breaking factor (used in Ru update)

%% Basic checks (lightweight)
if any([Htoe, Tm] <= 0) | any([Rc, htoe] < 0) | any([alpha, m] <= 0)
    error('Inputs must satisfy: Htoe>0, Tm>0, alpha>0, m>0, and Rc, htoe >= 0.');
end

%% Helper: degrees form of slopes for tand/atand consistency 
alpha_deg = atand(alpha);  % structure slope angle in degrees
m_deg     = atand(m);      % foreshore slope angle in degrees

%% Deep-water wavelength and Iribarren skeleton
L0 = 1.56 .* Tm.^2;        % [m] deep-water wavelength (linear theory approx)

%% Iteration: Equivalent slope via run-up Ru (Altomare 2016)
Ru      = 1.5 .* Htoe;    % initial guess for run-up (m)
tol_pc  = 1.0;             % [%] relative change tolerance
max_it  = 100;             % safety cap on iterations
it      = 0;

while true
    % Composite horizontal length (Eq. 13 in Altomare et al. 2016)
    Lslope = (Ru + htoe) ./ tand(alpha_deg) + (1.5*Htoe - htoe) ./ tand(m_deg);

    % Equivalent slope (Eq. 11)
    tan_delta = (1.5.*Htoe + Ru) ./ Lslope;

    % Iribarren with equivalent slope
    Ir_delta = tan_delta ./ sqrt(Htoe ./ L0);

    % Updated run-up (Eq. 12 with capping)
    Ru_new = min( 1.65 .* gamma_break .* Htoe .* Ir_delta, ...
              (4.0 - 1.5 ./ sqrt(Ir_delta)) .* gamma_break .* Htoe );  %% CHANGE

    it = it + 1;
    if it > max_it
        warning('q_altomare:NoConvergence', ...
            'Run-up iteration did not reach 1%% tolerance within %d steps.', max_it);
        Ru = Ru_new;
        break
    end

    % Convergence check (use absolute relative change)
    rel_change_pc = abs((Ru_new - Ru) ./ Ru) * 100;
    Ru = Ru_new;
    if rel_change_pc <= tol_pc
        break
    end
end

%% Mean discharge: dimensionless Q and dimensional q 
% Altomare regression intercept and scatter
c        = -0.791;     % intercept in base-10 logarithmic regression
sigma_c  = 0.294;      % 1-sigma scatter of intercept
z_95ones = 1.64;       % ~95% one-sided normal quantile

% Dimensionless mean discharge (Q)
denom = Htoe .* (0.33 + 0.022 .* Ir_delta);       % Htoe * scale(ξ_δ)
Q     = (10.^c) .* exp( - Rc ./ denom );

fac95 = 10.^( z_95ones * sigma_c);   % ≈ 3.0350
fac5  = 10.^(-z_95ones * sigma_c);   % ≈ 0.3295

Q_95 = fac95 .* Q;
Q_5  = fac5  .* Q;

q_95 = Q_95 .* sqrt(g .* Htoe.^3);
q_5  = Q_5  .* sqrt(g .* Htoe.^3);
q = Q .* sqrt(g .* Htoe.^3);

end




