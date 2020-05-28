clear all
close all
clc

load('air_properties_locked_data.mat','rho_0', 'v_sound', 'visc');
%% Parametric Sweep Space
pressure_sweep = logspace(log10(1000), log10(25000), 30);      % Pressure amplitude in Pa
frequency_sweep = [20e03];    % Acoustic Frequency in Hz
update_freq_sweep = [10];
r_min = 0.1;                      % min normalised radius
r_max = 0.2;                          % max nomalised radius
r_incr = 400;                       % number of increments in radius
rho_min = rho_0;                     % min density in kg m3
rho_max = 15;                     % max density in kg m3
rho_incr = 400;                     % number of increments in density
save_data = 0;
r_normalised = linspace(r_min, r_max, r_incr);
rho_screen = linspace(rho_min, rho_max, rho_incr);
[RR, RHORHO]= ndgrid(r_normalised, rho_screen);
acoustic_f = frequency_sweep;
%% Export Data Matrix
r_r = zeros(1, length(pressure_sweep));
rho_r = zeros(1, length(pressure_sweep));
max_scr = zeros(1, length(pressure_sweep));
vel_r = zeros(1, length(pressure_sweep));
omega_r = zeros(1, length(pressure_sweep));
tc_interest = [18 20 22]; % interested index
plot_count = 1;

for ii = 1:length(pressure_sweep)
        if sum(tc_interest == ii) == 1
        p_a = pressure_sweep(ii);
        [omega_n_store, max_velocity_store, drag_store] = main_function_fig5(p_a,acoustic_f, rho_0, v_sound, visc, r_min, r_max, r_incr, rho_min, rho_max, rho_incr, save_data);
        
        f_r = update_freq_sweep;
        wavelength = v_sound/acoustic_f;
        
        %% Sec 1 Screen Size
        Sc_size_sec1 = 2.*sqrt((max_velocity_store'.*(RR.*wavelength))./(pi.*f_r));
        omega_sec_1 = sqrt((f_r.*pi.*max_velocity_store')./(RR.*(wavelength)));
        Sc_size_sec1(omega_sec_1>omega_n_store') = NaN;
        
        %% Sec 2 Screen Size
        vel_sec2 = (((( max_velocity_store'.*omega_n_store').^2).*RR.*wavelength)./(f_r.*pi)).^(1/3);
        vel_sec2(vel_sec2 > max_velocity_store') = NaN;
        Sc_size_sec2 = 2.*sqrt((vel_sec2.*(RR.*wavelength))./(pi.*f_r));
        omega_sec_2 = (1./(vel_sec2).*max_velocity_store'.*omega_n_store');
        
        %% Get maxima of either sec
        sc_max_sec1 = max(max(Sc_size_sec1));
        sc_max_sec2 = max(max(Sc_size_sec2));
        
        %% Export
        if sc_max_sec1 > sc_max_sec2
            max_scr(ii) = sc_max_sec1;
            find_p = find(Sc_size_sec1==sc_max_sec1);
            mx_o = omega_sec_1(find_p);
            tt = max_velocity_store';
            mx_v = tt(find_p);
        else
            max_scr(ii) = sc_max_sec2;
            find_p = find(Sc_size_sec2==sc_max_sec2);
            mx_o = omega_sec_2(find_p);
            mx_v = vel_sec2(find_p);
        end
        
        mx_r = RR(find_p);
        mx_rho = RHORHO(find_p);
        
        r_r(ii) = mx_r;
        rho_r(ii) = mx_rho;
        vel_r(ii) = mx_v;
        omega_r(ii) = mx_o;
        
        %
%         if sum(tc_interest == ii) == 1
        figure(1)
        subplot(3,1,plot_count)
        hold on
        try %needed because there are points where there are no damping limited region
            [C, h] = contourf(rho_screen./rho_0, r_normalised, Sc_size_sec1.*1000,[0:2.5:100]);
        catch
            
        end
        [Cv, hv] = contourf(rho_screen./rho_0, r_normalised, Sc_size_sec2.*1000,[0:2.5:100]);
        scatter(mx_rho./rho_0, mx_r,60, 'rx','LineWidth',3)
        title([num2str(p_a) ' Pa'])
        ylabel('r/\lambda [-]')
        xlabel('\rho_p/\rho_0 [-]')
        set(gca,'FontSize',20)
        
        plot_count = plot_count + 1;        
    end
end
h123 = colorbar('Position', [0.6916+0.2134+0.02  0.1100+0.09  0.025  0.7], ...
    'YTick',[0:2.5:100],'FontSize',30);
title(h123, 'S_r [mm]','FontSize',25);
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'results\supp_figure\liftoff_1.png')

figure('units','normalized','outerposition',[0 0 1 1])
ha = tight_subplot(3,1,[.04 .04],[.15 .01],[.1 .05]);

%% Screen Size
axes(ha(1))
hold on
plot(pressure_sweep, max_scr.*1e03,'-ok','LineWidth',3)
scatter(pressure_sweep(tc_interest), max_scr(tc_interest).*1e03,200,'ro','LineWidth',3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'xticklabel',[])
yticks('auto')
yticklabels('auto')
grid on
grid minor
ylabel('S_r [mm]')
h = legend('20 kHz','Location','northwest');
h.FontSize = 25;
set(gca,'FontSize', 30)

%% Radius
axes(ha(2))
hold on
plot(pressure_sweep, r_r,'-ok','LineWidth',3)
scatter(pressure_sweep(tc_interest), r_r(tc_interest),200,'ro','LineWidth',3)
set(gca, 'XScale', 'log')
ylabel('r/\lambda [-]')
set(gca,'xticklabel',[])
yticks('auto')
yticklabels('auto')
grid on
grid minor
set(gca,'FontSize', 30)

%% Density
axes(ha(3))
hold on
plot(pressure_sweep, rho_r./rho_0,'-ok','LineWidth',3)
scatter(pressure_sweep(tc_interest), rho_r(tc_interest)./rho_0,200,'ro','LineWidth',3)
set(gca, 'XScale', 'log')
ylabel('\rho_p/\rho_0 [-]')
set(gca,'xticklabel',{'10^3', '10^4','10^5'})
yticks('auto')
yticklabels('auto')
grid on
grid minor
set(gca,'FontSize', 30)
xlabel('p_0 [Pa]','FontSize',30)

%% Export Figure
saveas(gcf,'results\supp_figure\liftoff_2.png')