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

function graph_gen_FIG5(calculate_parametric)
if calculate_parametric
    load('air_properties_locked_data.mat','rho_0', 'v_sound', 'visc');
    %% Parametric Sweep Space
    pressure_sweep = logspace(log10(1000), log10(25000), 30);      % Pressure amplitude in Pa
    frequency_sweep = [20e03 40e03 80e03 160e03];    % Acoustic Frequency in Hz
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
    
    %% Export Data Matrix
    r_r = zeros(length(frequency_sweep), length(pressure_sweep));
    rho_r = zeros(length(frequency_sweep), length(pressure_sweep));
    max_scr = zeros(length(frequency_sweep), length(pressure_sweep));
    vel_r = zeros(length(frequency_sweep), length(pressure_sweep));
    omega_r = zeros(length(frequency_sweep), length(pressure_sweep));
    tc_count = 1;
    for jj = 1:length(frequency_sweep)
        acoustic_f = frequency_sweep(jj);
        for ii = 1:length(pressure_sweep)
            p_a = pressure_sweep(ii);
            [omega_n_store, max_velocity_store] = main_function(p_a,acoustic_f, rho_0, v_sound, visc, r_min, r_max, r_incr, rho_min, rho_max, rho_incr, save_data);
            
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
                max_scr(jj, ii) = sc_max_sec1;
                find_p = find(Sc_size_sec1==sc_max_sec1);
                mx_o = omega_sec_1(find_p);
                tt = max_velocity_store';
                mx_v = tt(find_p);
            else
                max_scr(jj, ii) = sc_max_sec2;
                find_p = find(Sc_size_sec2==sc_max_sec2);             
                mx_o = omega_sec_2(find_p);
                mx_v = vel_sec2(find_p);
            end
            
            mx_r = RR(find_p);
            mx_rho = RHORHO(find_p);
            
            r_r(jj, ii) = mx_r;
            rho_r(jj, ii) = mx_rho;
            vel_r(jj, ii) = mx_v;
            omega_r(jj, ii) = mx_o;
        end
    end
    save('results\parametric_study.mat')
end
%% Images
load('results\parametric_study.mat') %#ok<LOAD>

figure('units','normalized','outerposition',[0 0 1 1])
ha = tight_subplot(3,1,[.04 .04],[.15 .01],[.1 .05]);

%% Screen Size
axes(ha(1))
hold on
plot(pressure_sweep, max_scr(1,:).*1e03,'-ok','LineWidth',3)
plot(pressure_sweep, max_scr(2,:).*1e03,'-^r','LineWidth',3)
plot(pressure_sweep, max_scr(3,:).*1e03,'--ob','LineWidth',3)
plot(pressure_sweep, max_scr(4,:).*1e03,':o','LineWidth',3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'xticklabel',[])
yticks('auto')
yticklabels('auto')
grid on
grid minor
ylabel('S_r [mm]')
h = legend('20 kHz', '40 kHz','80 kHz','160 kHz','Location','northwest');
h.FontSize = 25;
set(gca,'FontSize', 30)

%% Radius
axes(ha(2))
hold on
plot(pressure_sweep, r_r(1,:),'-ok','LineWidth',3)
plot(pressure_sweep, r_r(2,:),'-^r','LineWidth',3)
plot(pressure_sweep, r_r(3,:),'--ob','LineWidth',3)
plot(pressure_sweep, r_r(4,:),':o','LineWidth',3)
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
plot(pressure_sweep, rho_r(1,:)./rho_0,'-ok','LineWidth',3)
plot(pressure_sweep, rho_r(2,:)./rho_0,'-^r','LineWidth',3)
plot(pressure_sweep, rho_r(3,:)./rho_0,'--ob','LineWidth',3)
plot(pressure_sweep, rho_r(4,:)./rho_0,':o','LineWidth',3)
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
fig= gcf;
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[-0.125 -0.125 50 30],...
    'PaperSize',[50 30]);
print(fig,['results\main_figure\FIG5_APL_parametric'],'-dpdf','-r0')
end
