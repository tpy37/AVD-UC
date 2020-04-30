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

function [out] = spherical_Bessel_second_kind(nu,Z)
% Abramowitz and Stegun.
% Handbook of Mathematical Functions.
% page 437

out = sqrt((pi/2)./(Z)) .* bessely(nu+0.5, Z);

end

