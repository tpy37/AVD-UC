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

load('air_properties_locked_data.mat')

rho_p = 19;             % Density in kg m^{-3}
c_c = 900;              % Compressional wave speed in m s^{-1}
c_s =  900;             % Shear wave speed in m s^{-1}

acoustic_f = 40e03;     % Acoustic frequency in Hz
p_a = 5e03;             % Pressure amplitude in Pa

lambda = (v_sound/acoustic_f);
r = linspace(lambda/100, lambda,200);

figure('units','normalized','outerposition',[0 0 1 1])
hold on 
m_sum_tol = [1e-10];
F_rigid = abs(alpha_c(m_sum_calculator(r, acoustic_f, v_sound, rho_0, rho_p,m_sum_tol),p_a,rho_0, acoustic_f));
F_general = abs(alpha_c(m_sum_calculator_general(r, acoustic_f, v_sound, rho_0, rho_p, c_c, c_s, m_sum_tol),p_a,rho_0, acoustic_f));
plot(r./lambda, F_rigid,'o-','LineWidth',3)
plot(r./lambda, F_general,'x-','LineWidth',3)   

plot([1e-03 1e-03]./lambda, [0 6e-04],'r','LineWidth',3)

legend('Chen & Apfel (Rigid)', 'Chen & Apfel (c_c = 900, c_s = 900)','Point A','Location','southeast')
xlabel('r/\lambda [-]')
ylabel('Force [N]')
xlim([lambda/100 0.25])

% (interpn(r, F_rigid, 1e-03,'spline')./interpn(r, F_general, 1e-03,'spline')).*100

grid on
grid minor
set(gca,'FontSize', 24)

saveas(gcf,'results\supp_figure\solid_elastic.png')
