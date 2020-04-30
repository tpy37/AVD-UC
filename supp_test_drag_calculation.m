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


%% Compares 
Re = logspace(log10(0.01), log10(1e06), 2000);
figure('units','normalized','outerposition',[0 0 1 1])
C_d_z_flemmer = (24.*10.^(0.261.*Re.^(0.369)-0.105.*Re.^(0.431)-(0.124)./(1+(log10(Re)).^2)))./(Re);
C_d_z_flemmer(Re>3e05) = NaN;

plot(Re, C_d_z_flemmer, 'LineWidth', 3)
hold on
plot(Re, C_d(Re), 'LineWidth', 3)

plot([69914.5632 69914.5632], [0.001 1e09],'r', 'LineWidth', 3)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'FontSize',24)
xlabel('Re [-]')
ylabel('C_d [-]')
legend('Flemmer & Banks','Morrison', 'Max Reynolds Number Considered','Location','northwest')
ylim([0.01 1e09])
grid on
grid minor
saveas(gcf,'results\supp_figure\drag_comparison.png')