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

function [Dm] = calculate_D_general(k,r, rho_0, rho_p, kc, ks, m)
%% Chen, X. & Apfel, R. E. Radiation force on a spherical object in an axisymmetric wave field and its application to the calibration of high-frequency transducers. J. Acoust. Soc. Am. 99, 713 (1996).
%% Equation 4a and 5 in Chen & Apfel
ka = k.*r;
kc_a = kc.*r;
ks_a = ks.*r;

A_m_k_sa = (spherical_Bessel_first_kind(m, ks_a)-ks_a.*derivative_first_kind_spherical_bessel(m, ks_a))...
    ./(m.^2+m-2-((ks_a.^2)./2));
A_m_k_ca = (spherical_Bessel_first_kind(m, kc_a)-kc_a.*derivative_first_kind_spherical_bessel(m, kc_a))...
    ./(m.^2+m-2-((kc_a.^2)./2));

A = (rho_0.*(ks_a.^2))./(2.*rho_p);
B = 1./(m.^2+m-2-((ks_a.^2)./2));
C = ((((m.^2+m).*spherical_Bessel_first_kind(m, ks_a))./(spherical_Bessel_first_kind(m, ks_a)+A_m_k_sa))...
    +((kc_a.*derivative_first_kind_spherical_bessel(m, kc_a))./(A_m_k_ca)))...
    ./((((m.^2+m).*A_m_k_sa)./(spherical_Bessel_first_kind(m, ks_a)+2.*A_m_k_sa))-...
    ((spherical_Bessel_first_kind(m, kc_a)+2.*A_m_k_ca)./(A_m_k_ca)));

PHI_M = A.*B.*C;

Dm = (spherical_Bessel_first_kind(m, ka).*PHI_M - ka.*derivative_first_kind_spherical_bessel(m, ka))...
    ./(spherical_Bessel_second_kind(m, ka).*PHI_M - ka.*derivative_second_kind_spherical_bessel(m, ka));

end

