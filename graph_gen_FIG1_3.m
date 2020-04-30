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

close all

p_a = 5000;
acoustic_f = 40e03;

load(['results/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa.mat'])
[RR, RHORHO]= ndgrid(r_normalised, rho_screen);

%% Figure 1
figure('units','normalized','outerposition',[0 0 1 1])
[C, hvv] = contourf(rho_screen./rho_0, r_normalised, arf_store'.*1e03, [0:0.25:2]);
colormap hot
hold on
h = colorbar;
shading interp
set(gca, 'XScale', 'log')
xlabel('\rho_p / \rho_0 [-]')
ylabel('r/\lambda [-]')
caxis([0 2])
clabel(C, hvv, [0:0.25:2],'FontSize',20,'LabelSpacing',80000,'Color','white')
[model] = sub_fig1_code();
plot(model.sect1_rho, model.sect1_rr, 'g', 'LineWidth',3);
plot(model.sect2_rho, model.sect2_rr, 'g', 'LineWidth',3);
plot(model.sect3_rho, model.sect3_rr, 'g', 'LineWidth',3);
plot(model.sect4_rho, model.sect4_rr, 'g', 'LineWidth',3);
plot(model.sect5_rho, model.sect5_rr, 'g', 'LineWidth',3);
h1 = scatter(19./rho_0, 1e-03/(v_sound/acoustic_f), 'wo','LineWidth', 8);
text(19./rho_0+2.5, 1e-03/(v_sound/acoustic_f),'Point A','FontSize', 30,'Color','w')
set(gca,'FontSize', 40)

ylim([min(r_normalised) max(r_normalised)])

grid on
grid minor
title(h, 'mN')
set(gcf,'Units','inches');
set(gcf,...
    'PaperPosition',[-0.125 -0.125 50 30],...
    'PaperSize',[50 30]);
fig= gcf;
set(gca,'layer','top')
print(fig,['results\main_figure\FIG1_APL_arf_apfel' ],'-dpdf','-r0')
% return
%% Figure 2

figure('units','normalized','outerposition',[0 0 1 1])
[C, hvv] = contourf(rho_screen./rho_0, r_normalised, max_velocity_store',[0:2.5:20]);
shading interp
set(gca, 'XScale', 'log')
xlabel('\rho_p / \rho_0 [-]')
ylabel('r/\lambda [-]')
set(gca,'FontSize', 40)
grid on
grid minor
h = colorbar;
clabel(C, hvv, 0:5:20,'FontSize',25,'LabelSpacing',800,'Color','white')
caxis([0 20])
title(h, 'm s^{-1}')
hold on
h1 = scatter(19./rho_0, 1e-03/(v_sound/acoustic_f), 'ko','LineWidth', 8);

% Extract Points at Hirayama, R., Martinez Plasencia, D., Masuda, N. &
% Subramanian, S. A volumetric display for visual, tactile and audio 
% presentation using acoustic trapping. Nature 575, 320–323 (2019).

point_a_arf = interpn(rho_screen, r, arf_store, 19, 1e-03);
point_a_vel = interpn(rho_screen, r, max_velocity_store, 19, 1e-03);
% Output in Command
disp(['ARF at Point A = ' num2str(point_a_arf)]) 
disp(['VEL at Point A = ' num2str(point_a_vel)])

% Extract points for the maximum achievable velocity. 
best_p = max(max(max_velocity_store));
best_p_idx = find(max_velocity_store' == best_p);
best_p_r = RR(best_p_idx);
best_p_rho = RHORHO(best_p_idx)./rho_0;
% Output in Command
disp(['Best Velocity = ' num2str(best_p)])
disp(['r at Point O = ' num2str(best_p_r)])
disp(['rho at Point O = ' num2str(best_p_rho)])

% Draw Point A
text(19./rho_0+2.5, 1e-03/(v_sound/acoustic_f),'Point A','FontSize', 30)

% Draw Zero ARF Line
plot(model.sect1_rho, model.sect1_rr, 'g', 'LineWidth',3);
plot(model.sect2_rho, model.sect2_rr, 'g', 'LineWidth',3);
plot(model.sect3_rho, model.sect3_rr, 'g', 'LineWidth',3);
plot(model.sect4_rho, model.sect4_rr, 'g', 'LineWidth',3);
plot(model.sect5_rho, model.sect5_rr, 'g', 'LineWidth',3);

ylim([min(r_normalised) max(r_normalised)])

set(gcf,'Units','inches');
set(gcf,...
    'PaperPosition',[-0.125 -0.125 50 30],...
    'PaperSize',[50 30]);
fig= gcf;
set(gca,'layer','top')
% % Export as PDF
print(fig,['results\main_figure\FIG2_APL_max_v' ],'-dpdf','-r0')


%% Figure 3
figure('units','normalized','outerposition',[0 0 1 1])

pcolor(rho_screen./rho_0, r_normalised, log10(omega_n_store)')
shading interp
set(gca, 'XScale', 'log')
xlabel('\rho_p / \rho_0 [-]')
ylabel('r/\lambda [-]')
caxis([min(min(log10(omega_n_store))) max(max(log10(omega_n_store)))])
cb = colorbar('YTick',[-1:1:6],'YTickLabel',{'10^{-1}','10^0','10^1','10^2','10^3','10^4','10^5','10^6'});
set(gca,'FontSize', 40)
title(cb, 'rad s^{-1}')
grid on
grid minor
set(gca,'layer','top')
hold on

w_g = sqrt((f_r.*pi.*max_velocity_store')./(RR.*(v_sound/acoustic_f)));
idx = find(abs(w_g-omega_n_store')< 1.5);
w_g(w_g > omega_n_store') = NaN;
Q = bwperim(bwmorph(~isnan( w_g ),'thicken', 1)); % Using the contrast between the boundary to find the outline.
Q(1,:) = 0;
Q(end, :) = 0;
Q(:,1) = 0;
Q(:,end) = 0;
limit_points = find(Q == 1);
X = RHORHO(limit_points)./rho_0;
Y = RR(limit_points);

% As line was discontinuous, add more points to specific region
temp_r_min = 1/100; %min normalised radius
temp_r_max = 0.18; %max nomalised radius
temp_r_incr = 600; % number of increments in radius
temp_rho_min = min(rho_screen./rho_0); %min density in kg m3
temp_rho_max = 1; %max density in kg m3
temp_rho_incr = 400; % number of increments in density

[temp_omega, temp_max_v] = main_function(p_a,acoustic_f, rho_0, v_sound, visc, temp_r_min, temp_r_max, temp_r_incr, temp_rho_min, temp_rho_max, temp_rho_incr, 0);
[temp_rr, tmep_rhorho] = ndgrid(linspace(temp_r_min, temp_r_max, temp_r_incr), logspace(log10(temp_rho_min), log10(temp_rho_max), temp_rho_incr));
temp_w_g = sqrt((f_r.*pi.*temp_max_v')./(temp_rr.*(v_sound/acoustic_f)));

temp_w_g(temp_w_g > temp_omega') = NaN;
Q = bwperim(bwmorph(~isnan( temp_w_g ),'thicken', 1));
Q(1,:) = 0;
Q(end, :) = 0;
Q(:,1) = 0;
Q(:,end) = 0;
limit_points = find(Q == 1);
X_temp = tmep_rhorho(limit_points)./rho_0;
Y_temp = temp_rr(limit_points);

scatter_p = scatter([X; X_temp], [Y;Y_temp],'ro','LineWidth',2.5);

h1 = scatter(19./rho_0, 1e-03/(v_sound/acoustic_f), 'ko','LineWidth', 8);

text(19./rho_0+2.5, 1e-03/(v_sound/acoustic_f),'Point A','FontSize', 30);
% return
text(0.09, 0.3, 'v = v_{max}, a \leq a_{max}','FontSize', 34)
text(0.937, 0.54, 'v \leq v_{max}, a = a_{max}','FontSize', 34)

plot(model.sect1_rho, model.sect1_rr, 'g', 'LineWidth',3);
plot(model.sect2_rho, model.sect2_rr, 'g', 'LineWidth',3);
plot(model.sect3_rho, model.sect3_rr, 'g', 'LineWidth',3);
plot(model.sect4_rho, model.sect4_rr, 'g', 'LineWidth',3);
plot(model.sect5_rho, model.sect5_rr, 'g', 'LineWidth',3);
legend(scatter_p, 'Transition Point (10 Hz)')
ylim([min(r_normalised) max(r_normalised)])
shading interp
set(gcf,'Units','inches');
set(gcf,...
    'PaperPosition',[-0.125 -0.125 50 30],...
    'PaperSize',[50 30]);
fig= gcf;
set(gca,'layer','top')
print(fig,['results\main_figure\FIG3_omega_g_max' ],'-dpdf','-r0')

