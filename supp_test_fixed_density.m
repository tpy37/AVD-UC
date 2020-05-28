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

function supp_test_fixed_density(calculate_parametric)
if calculate_parametric
    load('air_properties_locked_data.mat','rho_0', 'v_sound', 'visc');
    %% Parametric Sweep Space
    pressure_sweep = logspace(log10(1000), log10(25000), 30);      % Pressure amplitude in Pa
    frequency_sweep = logspace(log10(20e03), log10(160e03), 40);    % Acoustic Frequency in Hz
    update_freq_sweep = [10];
    
    density_point_A = 19;
    r_min = 1/100;   % min normalised radius
    r_max = 1;                          % max nomalised radius
    r_incr = 600;                       % number of increments in radius
    rho_min = density_point_A;          % min density in kg m3
    rho_max = density_point_A;          % max density in kg m3
    rho_incr = 1;                       % number of increments in density
    save_data = 0;
    r_normalised = linspace(r_min, r_max, r_incr);
    rho_screen = logspace(log10(rho_min), log10(rho_max), rho_incr);
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
            wavelength = v_sound/acoustic_f;

            [omega_n_store, max_velocity_store] = main_function(p_a,acoustic_f, rho_0, v_sound, visc, r_min, r_max, r_incr, rho_min, rho_max, rho_incr, save_data);
            f_r = update_freq_sweep;
            
            
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
    save('results\supp_figure\sub_data\point_A_study.mat')
end
%% Images
load('results\supp_figure\sub_data\point_A_study.mat') %#ok<LOAD>

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
hold on
surf(frequency_sweep./1e03,pressure_sweep./1e03, max_scr'*1e03)
colormap parula
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-41.6500, 24.4808)
xlabel('Frequency [kHz]')
ylabel('p_0 [kPa]')
zlabel('S_r [mm]')
set(gca,'FontSize',24)
grid on
grid minor
subplot(1,2,2)
hold on
surf(frequency_sweep./1e03,pressure_sweep./1e03, r_r')
colormap parula
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
view(-41.6500, 24.4808)
xlabel('Frequency [kHz]')
ylabel('p_0 [kPa]')
zlabel('r/\lambda [-]')
set(gca,'FontSize',24)
grid on
grid minor
saveas(gcf,'results\supp_figure\fixed_density_plot.png')

end