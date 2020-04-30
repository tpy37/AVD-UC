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

% The main_function.m in this folder is programmed for rigid-assumption. In
% order to evaluate scattering from a specific material, you have to
% replace the "m_sum_calculator" function with "m_sum_calculator_general".
% The m_sum_calculator_general takes additional inputs, c_c and c_s which
% are compressional and shear wave speed, respectively. Some examples of
% using m_sum_calculator_general are shown below. alpha.m evaluates the
% amplitude of acoustic radiation force. 

clear all
close all
clc

load('air_properties_locked_data.mat')

rho_p = 10;             % Density in kg m^{-3}
c_c = 400;              % Compressional wave speed in m s^{-1}
c_s =  500;             % Shear wave speed in m s^{-1}

acoustic_f = 40e03;     % Acoustic frequency in Hz
p_a = 1;                % Pressure amplitude in Pa

lambda = (v_sound/acoustic_f);
r = linspace(lambda/100, lambda,200);

figure('units','normalized','outerposition',[0 0 1 1])
hold on 
m_sum_tol = [1e-10];
plot(r./lambda, abs(alpha(m_sum_calculator(r, acoustic_f, v_sound, rho_0, rho_p,m_sum_tol),p_a,rho_0, acoustic_f)),'o','LineWidth',3)
plot(r./lambda, abs(alpha(m_sum_calculator_general(r, acoustic_f, v_sound, rho_0, rho_p, c_c, c_s, m_sum_tol),p_a,rho_0, acoustic_f)),'LineWidth',3)   

legend('Chen & Apfel (Rigid)', 'Chen & Apfel (c_c = 400, c_s = 500)')
xlabel('r/\lambda [-]')
ylabel('Force [N]')
ylim([-1e-10 5e-10])
grid on
grid minor
set(gca,'FontSize', 24)