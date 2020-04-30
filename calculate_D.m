% Acoustophoretic Volumetric Display - Ultimate Capability (AVD-UC).
% 
% Codes for: "What is the Ultimate Capability of Acoustophoretic Volumetric
% Displays?".
% 
% Authors: Tatsuki Fushimi, Bruce W. Drinkwater, and Thomas L. Hill.
% 
% Journal: Applied Physics Letters.
% 
% Created: 30-Apr-2020.
% 
% Please contact Tatsuki FUSHIMI (t.fushimi@bristol.ac.uk) for any inquiry.
% Released under MIT License

function [Dm] = calculate_D(k,r, rho_0, rho_p, m)
%% Chen, X. & Apfel, R. E. Radiation force on a spherical object in an axisymmetric wave field and its application to the calibration of high-frequency transducers. J. Acoust. Soc. Am. 99, 713 (1996).
%% Equation 4c in Chen & Apfel
ka = k.*r;
if m == 1
    rho_ratio = rho_0/rho_p;
    Dm = (rho_ratio*spherical_Bessel_first_kind(m, ka)-ka.*derivative_first_kind_spherical_bessel(m, ka))./(rho_ratio.*spherical_Bessel_second_kind(m, ka)-ka.*derivative_second_kind_spherical_bessel(m, ka));
else
    Dm = derivative_first_kind_spherical_bessel(m, ka)./derivative_second_kind_spherical_bessel(m, ka);
end
end

