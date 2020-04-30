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

clear all
close all
clc

load('air_properties_locked_data.mat','rho_0', 'v_sound', 'visc');

p_a = [5e03];      % Pressure amplitude in Pa
acoustic_f = [40e03];    % Acoustic Frequency in Hz
r_min = 1/100;                      % min normalised radius
r_max = 1;                          % max nomalised radius
r_incr = 600;                       % number of increments in radius
rho_min = 1e-1;                     % min density in kg m3
rho_max = 1000;                     % max density in kg m3
rho_incr = 400;                     % number of increments in density
save_data = 0;                      % Export data to .mat file

%% Run Analysis
[~, max_velocity_store_ver] = main_function(p_a,acoustic_f, rho_0, v_sound, visc, r_min, r_max, r_incr, rho_min, rho_max, rho_incr, save_data);
[~, max_velocity_store_hoz] = main_function_horizontal(p_a,acoustic_f, rho_0, v_sound, visc, r_min, r_max, r_incr, rho_min, rho_max, rho_incr, save_data);

r_normalised = linspace(r_min, r_max, r_incr);
rho_screen = logspace(log10(rho_min), log10(rho_max), rho_incr);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
[C, hvv] = contourf(rho_screen./rho_0, r_normalised, max_velocity_store_ver',[0:2.5:20]);
shading interp
set(gca, 'XScale', 'log')
xlabel('\rho_p / \rho_0 [-]')
ylabel('r/\lambda [-]')
set(gca,'FontSize', 20)
grid on
grid minor
h = colorbar;
clabel(C, hvv, 0:5:20,'FontSize',25,'LabelSpacing',800,'Color','white')
caxis([0 20])
title(h, 'm s^{-1}')
hold on
h1 = scatter(19./rho_0, 1e-03/(v_sound/acoustic_f), 'ko','LineWidth', 8);

% Draw Point A
text(19./rho_0+2.5, 1e-03/(v_sound/acoustic_f),'Point A','FontSize', 30)
ylim([min(r_normalised) max(r_normalised)])
title('Vertical (W/ grav/buoy)')
subplot(1,2,2)
[C, hvv] = contourf(rho_screen./rho_0, r_normalised, max_velocity_store_hoz',[0:2.5:20]);
shading interp
set(gca, 'XScale', 'log')
xlabel('\rho_p / \rho_0 [-]')
ylabel('r/\lambda [-]')
set(gca,'FontSize', 20)
grid on
grid minor
h = colorbar;
clabel(C, hvv, 0:5:20,'FontSize',25,'LabelSpacing',800,'Color','white')
caxis([0 20])
title(h, 'm s^{-1}')
hold on
h1 = scatter(19./rho_0, 1e-03/(v_sound/acoustic_f), 'ko','LineWidth', 8);

% Draw Point A
text(19./rho_0+2.5, 1e-03/(v_sound/acoustic_f),'Point A','FontSize', 30)
ylim([min(r_normalised) max(r_normalised)])
title('Horizontal (No grav/buoy)')
saveas(gcf,'results\supp_figure\vertical_horizontal_compare.png')