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

for ii = 1:length(pressure_sweep)
    p_a = pressure_sweep(ii);
    for jj = 1:length(frequency_sweep)
        acoustic_f = frequency_sweep(jj);
        main_function(p_a,acoustic_f, rho_0, v_sound, visc, r_min, r_max, r_incr, rho_min, rho_max, rho_incr, save_data);
    end
end