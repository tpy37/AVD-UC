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

rho_p = 40;
acoustic_f = 40e03;
p_a = 1;
lambda = (v_sound/acoustic_f);
r = linspace(lambda/100, lambda,200);

figure('units','normalized','outerposition',[0 0 1 1])
hold on
%% Gorkov (Rigid assumption)
E_gorkov = (p_a.^2)./(4.*rho_0*(v_sound.^2));
f_1 = 1;
f_2 = 1;
PHI = 1/3*f_1 + 0.5*f_2; 
F_gorkov = 4.*pi.*PHI.*(2*pi*acoustic_f./v_sound).*(r.^3).*E_gorkov;
subplot(1,2,1)
hold on
plot(r./lambda, F_gorkov, 'LineWidth',3)

%% Chen & Apfel
m_sum_tol_screen = [1e-10 1e-30];

l_cor = {'+','d'};

for ii = 1:length(m_sum_tol_screen)
    m_sum_tol = m_sum_tol_screen(ii);
    plot(r./lambda, abs(alpha_c(m_sum_calculator(r, acoustic_f, v_sound, rho_0, rho_p, m_sum_tol),p_a,rho_0, acoustic_f)),l_cor{ii},'LineWidth',1)
end
    
legend('Gorkov', 'Chen & Apfel (1\times10^{-10})', 'Chen & Apfel (1\times10^{-30})')
xlabel('r/\lambda [-]')
ylabel('Force [N]')
ylim([-1e-10 2e-10])
grid on
grid minor
set(gca,'FontSize', 24)

subplot(1,2,2)
hold on 
plot(r./lambda, F_gorkov, 'LineWidth',3)
m_sum_tol_screen = [1e-10 1e-30];

l_cor = {'+','d'};

for ii = 1:length(m_sum_tol_screen)
    m_sum_tol = m_sum_tol_screen(ii);
    plot(r./lambda, abs(alpha_c(m_sum_calculator(r, acoustic_f, v_sound, rho_0, rho_p, m_sum_tol),p_a,rho_0, acoustic_f)),l_cor{ii},'LineWidth',2)
end
    
legend('Gorkov', 'Chen & Apfel (1\times10^{-10})', 'Chen & Apfel (1\times10^{-30})')
xlabel('r/\lambda [-]')
ylabel('Force [N]')
ylim([-1e-10 1e-10])
xlim([0 0.2])
grid on
grid minor
set(gca,'FontSize', 24)

saveas(gcf,'results\supp_figure\arf_comparison.png')