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

function supp_test_force_sepp(recalculate)
load('air_properties_locked_data.mat','rho_0', 'v_sound', 'visc');

pressure_sweep = [5e03 25e03];      % Pressure amplitude in Pa
frequency_sweep = [20e03 40e03];    % Acoustic Frequency in Hz
%% Set Rendering Frequency
f_r = 10;
if recalculate
    tc_count = 1;
    cyc_len = 300;
    timing_function_raphson = @(w_p_i) linspace(0, 1./(w_p_i./(2*pi)), cyc_len);
    vel_fun = @(v_amp, w_p_i) - v_amp.*sin(w_p_i.*timing_function_raphson(w_p_i));
    

    S_r_ratio = linspace(0.8, 1.05, 20);
    for ii = 1:length(pressure_sweep)
        p_a = pressure_sweep(ii);
        for jj = 1:length(frequency_sweep)
            acoustic_f = frequency_sweep(jj);
            load(['results/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa.mat'])
            
            if tc_count == 1
                [RR, RHORHO]= ndgrid(r_normalised, rho_screen);
            end
            
            wavelength = v_sound/acoustic_f;
            Sc_size_sec1 = 2.*sqrt((max_velocity_store'.*(RR.*wavelength))./(pi.*f_r));
            
            omega_g_v = sqrt((f_r.*pi.*max_velocity_store')./(RR.*(wavelength)));
            Sc_size_sec1(omega_g_v>omega_n_store') = NaN;
                      
            v_sec2 = (((( max_velocity_store'.*omega_n_store' ).^2).*RR.*wavelength)./(f_r.*pi)).^(1/3);
            v_sec2(v_sec2 > max_velocity_store') = NaN;
            
            omega_sec2 = (1./(v_sec2).*max_velocity_store'.*omega_n_store');
            Sc_size_sec2 = 2.*sqrt((v_sec2.*(RR.*wavelength))./(pi.*f_r));
            ratio_store_sec1 = zeros(length(r_normalised), length(rho_screen));
            ratio_store_sec2 = zeros(length(r_normalised), length(rho_screen));
            
            permissible_SC_size = zeros(length(rho_screen), length(r_normalised));
            permissible_v = zeros(length(rho_screen), length(r_normalised));
            permissible_o = zeros(length(rho_screen), length(r_normalised));
            permissible_ratio = zeros(length(rho_screen), length(r_normalised));
            
            for kk = 1:length(rho_screen)
                rho_p = rho_screen(kk);
                
                for ll = 1:length(r_normalised)
                    
                    Mass_total = mass(r_normalised(ll).*wavelength,rho_p)+virtual_mass(r_normalised(ll).*wavelength,rho_0); % inertia needs to consider both particle and added mass
                    drag_fun = @(v_amp, w_p_i) Drag_Force_Zero_Escape(2, r_normalised(ll).*wavelength,rho_0,vel_fun(v_amp, w_p_i),visc);
                    inert_fun = @(v_amp, w_p_i) -(-v_amp.*w_p_i.*cos(w_p_i.*timing_function_raphson(w_p_i))).*Mass_total;
                    buoy_grav_fun = F_gravity_buoyancy(r_normalised(ll).*wavelength,rho_p,rho_0);
                    total_force_f = @(v_c, o_c) buoy_grav_fun.*ones(1,cyc_len)+inert_fun(v_c, o_c)+drag_fun(v_c, o_c);
                    arf_amp =arf_store(kk, ll);
                    cost = @(v_c, o_c) max(abs(total_force_f(v_c, o_c)./arf_amp));
                    
                    omega_n_max = omega_g_v(ll, kk);
                    
                    v_cand = max_velocity_store(kk, ll);
                    total_force_sec1 = total_force_f(v_cand, omega_n_max);
                    omega_cand = omega_sec2(ll, kk);
                    total_force_sec2 = total_force_f(v_sec2(ll,kk), omega_cand);
                    
                    ratio_store_sec1(ll, kk) = max(abs(total_force_sec1./arf_amp));
                    ratio_store_sec2(ll, kk) = max(abs(total_force_sec2./arf_amp));
                    
                    if isnan(omega_sec2(ll,kk))
                        sc_choice = Sc_size_sec1(ll,kk);
                    else
                        sc_choice = Sc_size_sec2(ll,kk);
                    end
                    
                    if ~isnan(sc_choice)
                        o_output = zeros(1, length(S_r_ratio));
                        v_output = zeros(1, length(S_r_ratio));
                        sr_output = zeros(1, length(S_r_ratio));
                        cost_output = zeros(1,length(S_r_ratio));
                        for TT = 1:length(S_r_ratio)
                            
                            A_cand = sc_choice*S_r_ratio(TT)/2;
                            omega_cand = f_r*pi*A_cand/(r_normalised(ll).*wavelength);
                            v_cand = omega_cand*A_cand;
                            
                            o_output(TT) = omega_cand;
                            v_output(TT) = v_cand;
                            sr_output(TT) = sc_choice;
                            cost_output(TT) = cost(v_cand, omega_cand);
                        end
                        permissible_SC_size(kk,ll) = interpn(cost_output, sr_output, 1,'spline');
                        permissible_v(kk,ll) = interpn(cost_output, v_output, 1,'spline');
                        permissible_o(kk,ll) = interpn(cost_output, o_output, 1,'spline');
                        permissible_ratio(kk,ll) = cost(interpn(cost_output, v_output, 1,'spline'), interpn(cost_output, o_output, 1,'spline'));
                    else
                        permissible_SC_size(kk,ll) = NaN;
                        permissible_v(kk,ll) = NaN;
                        permissible_o(kk,ll) = NaN;
                        permissible_ratio(kk,ll) = NaN;
                    end
                end
            end
            SS_size = 2.*sqrt((permissible_v'.*(RR.*wavelength))./(pi.*f_r));
            
            save(['results/supp_figure/sub_data/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa.mat'],...
                'omega_g_v','permissible_SC_size','permissible_v','permissible_o','SS_size','ratio_store_sec1','ratio_store_sec2','v_sec2','Sc_size_sec1','Sc_size_sec2','permissible_ratio','omega_sec2')
        end
    end
end

figure
tc_count = 1;
for ii = 1:length(pressure_sweep)
    p_a = pressure_sweep(ii);
    for jj = 1:length(frequency_sweep)
        acoustic_f = frequency_sweep(jj);
        load(['results/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa.mat'])
        load(['results/supp_figure/sub_data/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa.mat'])
        ratio_store_sec1(omega_g_v>omega_n_store') = NaN;
        if tc_count == 1
            [RR, RHORHO]= ndgrid(r_normalised, rho_screen);
        end
        
        %% Figure for Ratio Contour
        figure(1)
        subplot(2,2,tc_count)
        hold on
        pcolor(rho_screen./rho_0, r_normalised, ratio_store_sec1)
        shading interp
        caxis([1 1.3])
        pcolor(rho_screen./rho_0, r_normalised, ratio_store_sec2)
        caxis([1 1.3])
        shading interp
        set(gca, 'XScale', 'log')
        ylim([min(r_normalised) max(r_normalised)])
        xlim([min(rho_screen)./rho_0 max(rho_screen)./rho_0])
        cbh = colorbar ; %Create Colorbar
        cbh.Ticks = linspace(1, 1.3, 2) ; %Create 8 ticks from zero to 1
        cbh.TickLabels = {'\leq 1','1.3'} ;
        xlabel('\rho_p / \rho_0 [-]')
        ylabel('r/\lambda [-]')
        grid on
        grid minor
        set(gca,'layer','top')
        set(gca,'FontSize',20)
        if tc_count == 1
            title('5 kPa, 20 kHz')
        elseif tc_count == 2
            title('5 kPa, 40 kHz')
        elseif tc_count == 3
            title('25 kPa, 20 kHz')
        else
            title('25 kPa, 40 kHz')
        end
        %% Figure for recalculated screen size
        figure(3)
        subplot(2,2,tc_count)
        hold on
        contourf(rho_screen./rho_0, r_normalised, SS_size.*1000)
        xlabel('\rho_p / \rho_0 [-]')
        ylabel('r/\lambda [-]')
        shading interp
        colorbar
        set(gca, 'XScale', 'log')
        grid on
        grid minor
        set(gca,'layer','top')
        set(gca,'FontSize',20)
        if tc_count == 1
            title('5 kPa, 20 kHz')
        elseif tc_count == 2
            title('5 kPa, 40 kHz')
        elseif tc_count == 3
            title('25 kPa, 20 kHz')
        else
            title('25 kPa, 40 kHz')
        end
        
        %% Disp
        find_sec1 = find(ratio_store_sec1>1.15);
        find_sec2 = find(ratio_store_sec2>1.15);
        compare_ss1 = Sc_size_sec1./SS_size;
        compare_ss2 = Sc_size_sec2./SS_size;
        
        m1 = min(compare_ss1(find_sec1));
        m2 = min(compare_ss2(find_sec2));
        m3 = max(compare_ss1(find_sec1));
        m4 = max(compare_ss2(find_sec2));
        disp(['------- Subplot ' num2str(tc_count) ' -------'])
        disp(['minimum diff in 15% boundary = ' num2str(min([m1 m2]))]);
        disp(['maximum diff in 15% boundary = ' num2str(max([m3 m4]))]);
        
        if tc_count == 2
            %% Plot Cross Section of Interest
            figure(1)
            plot([min(rho_screen)./rho_0 max(rho_screen)./rho_0],[0.35 0.35],'r','LineWidth',3)
            %% Cross Section Plot
            figure(2)
            %% Ratio
            subplot(4,1,1)
            hold on
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, ratio_store_sec1', rho_screen, 0.35),'-k','LineWidth',2);
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, permissible_ratio, rho_screen, 0.35),'--r','LineWidth',2);
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, ratio_store_sec2', rho_screen, 0.35),'-k','LineWidth',2);
            set(gca, 'XScale', 'log')
            xlim([min(rho_screen)./rho_0 max(rho_screen)./rho_0])
            ylabel('R')
            legend('Force Separation Method (Manuscript)','Correction')
            set(gca,'FontSize',24)
            subplot(4,1,2)
            %% Velocity
            hold on
            tt = max_velocity_store;
            tt(omega_g_v'>omega_n_store) = NaN;
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, tt, rho_screen, 0.35),'-k','LineWidth',2);
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, v_sec2', rho_screen, 0.35),'-k','LineWidth',2);
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, permissible_v, rho_screen, 0.35),'--r','LineWidth',2);
            set(gca, 'XScale', 'log')
            xlim([min(rho_screen)./rho_0 max(rho_screen)./rho_0])
            ylabel('v [m s^{-1}]')
            set(gca,'FontSize',24)
            subplot(4,1,3)
            %% Omega
            hold on
            tt = omega_g_v';
            tt(omega_g_v'>omega_n_store) = NaN;
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, tt, rho_screen, 0.35),'-k','LineWidth',2);
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, omega_sec2', rho_screen, 0.35),'-k','LineWidth',2);
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, permissible_o, rho_screen, 0.35),'--r','LineWidth',2);
            set(gca, 'XScale', 'log')
            xlim([min(rho_screen)./rho_0 max(rho_screen)./rho_0])
            ylabel('\omega_p [rad s^{-1}]')
            set(gca,'FontSize',24)
            subplot(4,1,4)
            %% Screen Size
            hold on
            tt = Sc_size_sec1';
            tt(omega_g_v'>omega_n_store) = NaN;
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, tt, rho_screen, 0.35).*1000,'-k','LineWidth',2);
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, Sc_size_sec2', rho_screen, 0.35).*1000,'-k','LineWidth',2);
            plot(rho_screen./rho_0, interpn(rho_screen, r_normalised, SS_size', rho_screen, 0.35).*1000,'--r','LineWidth',2);
            set(gca, 'XScale', 'log')
            xlim([min(rho_screen)./rho_0 max(rho_screen)./rho_0])
            ylabel('S_r [mm]')
            xlabel('\rho_p/\rho_0 [-]')
            set(gca,'FontSize',24)
            
            %% Display diff
            c_d1 = max(interpn(rho_screen, r_normalised, Sc_size_sec1', rho_screen, 0.35)./interpn(rho_screen, r_normalised, SS_size', rho_screen, 0.35));
            c_d2 = max(interpn(rho_screen, r_normalised, Sc_size_sec2', rho_screen, 0.35)./interpn(rho_screen, r_normalised, SS_size', rho_screen, 0.35));
            
            disp(['Max Difference in Plot 2 = ' num2str(max([c_d1 c_d2]))])
            
        end
        tc_count = tc_count +1;
    end
end
%% Export Figure
figure(1)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'results\supp_figure\overshoot_graph.png')
figure(2)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'results\supp_figure\overshoot_graph_anlysis.png')
figure(3)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'results\supp_figure\new_screen_size.png')
end



