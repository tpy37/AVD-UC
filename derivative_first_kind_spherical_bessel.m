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

function [out] = derivative_first_kind_spherical_bessel(nu,Z)
%https://www.wolframalpha.com/input/?i=first+derivative+of+spherical+bessel+function+of+first+kind

out = 0.5.*(spherical_Bessel_first_kind(nu-1, Z)-...
    ((spherical_Bessel_first_kind(nu, Z)+Z.*spherical_Bessel_first_kind(nu+1, Z))./(Z)));
end

