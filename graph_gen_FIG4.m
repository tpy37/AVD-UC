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

tc_count = 1;

figure('units','normalized','outerposition',[0 0 1 1])
ha = tight_subplot(2,2,[.04 .04],[.15 .08],[.1 .12]);

v_h = [0:40:260];

cyc_len = 300;
%% Go through the sweep
for ii = 1:length(pressure_sweep)
    p_a = pressure_sweep(ii);
    for jj = 1:length(frequency_sweep)
        acoustic_f = frequency_sweep(jj);
        load(['results/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa.mat'])
        
        if tc_count == 1
            [RR, RHORHO]= ndgrid(r_normalised, rho_screen);
        end
        
        wavelength = v_sound/acoustic_f;
        
        %% Sec 1 Screen Size
        Sc_size_sec1 = 2.*sqrt((max_velocity_store'.*(RR.*wavelength))./(pi.*f_r));
        omega_sec1 = sqrt((f_r.*pi.*max_velocity_store')./(RR.*(wavelength)));
        Sc_size_sec1(omega_sec1>omega_n_store') = NaN;
        
        %% Sec 2 Screen Size
        v_sec2 = (((( max_velocity_store'.*omega_n_store').^2).*RR.*wavelength)./(f_r.*pi)).^(1/3);
        v_sec2(v_sec2 > max_velocity_store') = NaN;
        omega_ac = (1./(v_sec2).*max_velocity_store'.*omega_n_store');
        Sc_size_sec2 = 2.*sqrt((v_sec2.*(RR.*wavelength))./(pi.*f_r));          
        
        %% Calculate Overshoot
        ratio_store_sec1 = zeros(length(rho_screen), length(r_normalised));
        ratio_store_sec2 = zeros(length(rho_screen), length(r_normalised));
        for kk = 1:length(rho_screen)
            for ll = 1:length(r_normalised)
                rho_p = rho_screen(kk);
                Mass_total = mass(r_normalised(ll).*wavelength,rho_p)+virtual_mass(r_normalised(ll).*wavelength,rho_0); % inertia needs to consider both particle and added mass
                
                timing_function_raphson = @(w_p_i) linspace(0, 1./(w_p_i./(2*pi)), cyc_len);
                vel_fun = @(v_amp, w_p_i) - v_amp.*sin(w_p_i.*timing_function_raphson(w_p_i));
                drag_fun = @(v_amp, w_p_i) Drag_Force_Zero_Escape(2, r_normalised(ll).*wavelength,rho_0,vel_fun(v_amp, w_p_i),visc);
                inert_fun = @(v_amp, w_p_i) -(-v_amp.*w_p_i.*cos(w_p_i.*timing_function_raphson(w_p_i))).*Mass_total;
                buoy_grav_fun = F_gravity_buoyancy(r_normalised(ll).*wavelength,rho_p,rho_0);
                omega_n_max = omega_sec1(ll, kk);
                v_cand = max_velocity_store(kk, ll);
                total_force_sec1 = buoy_grav_fun.*ones(1,cyc_len)+inert_fun(v_cand, omega_n_max)+drag_fun(v_cand, omega_n_max);
                omega_cand = omega_ac(ll,kk);
                total_force_sec2 = buoy_grav_fun.*ones(1,cyc_len)+inert_fun(v_sec2(ll,kk), omega_cand)+drag_fun(v_sec2(ll,kk), omega_cand);
                arf_amp =arf_store(kk, ll);
                
                ratio_store_sec1(kk, ll) = max(abs(total_force_sec1./arf_amp));
                ratio_store_sec2(kk, ll) = max(abs(total_force_sec2./arf_amp));
            end
        end
        ratio_store_sec1(omega_sec1'>omega_n_store) = NaN;
        ratio_store_sec2(v_sec2>=max_velocity_store') = NaN;
        
        
        %% Plot Contour Plots
        axes(ha(tc_count))
        
        hold on
        [C, h] = contourf(rho_screen./rho_0, r_normalised, Sc_size_sec1.*1000,[0:10:260]);
        [Cv, hv] = contourf(rho_screen./rho_0, r_normalised, Sc_size_sec2.*1000,[0:10:260]);
        hold on
        shading interp
        set(gca, 'XScale', 'log')
        
        %% Information Output
        temp_SC_SIZE = max(max(Sc_size_sec1));
        temp_SC_SIZE_AC = max(max(Sc_size_sec2));
        max_A_r = max([temp_SC_SIZE temp_SC_SIZE_AC]);
        if temp_SC_SIZE > temp_SC_SIZE_AC
            find_p = find(Sc_size_sec1==temp_SC_SIZE);
        else
            find_p = find(Sc_size_sec2==temp_SC_SIZE_AC);
        end
        disp(['------- Subplot ' num2str(tc_count) ' -------'])
        disp(['Best SCR Size = ' num2str(max_A_r)])
        disp(['rho @ best = ' num2str(RHORHO(find_p)./rho_0)])
        disp(['r @ best = ' num2str(RR(find_p))])
        
        point_a_acc = interpn(rho_screen, r, max_velocity_store.*omega_n_store, 19, 1e-03);
        point_a_omega = interpn(rho_screen, r, omega_n_store, 19, 1e-03);
        
        disp(['Max Acceleration @ Point A= ' num2str(max(max(point_a_acc)))])
        disp(['Max Particle Frequency @ Point A= ' num2str(max(max(point_a_omega)))])
        
        
        rem_ratio = abs(rem_check./ arf_store).*100;
        disp(['Largest abs(f_val/ARF) is (in percentage)= ' num2str(max(max(rem_ratio(~isnan(rem_ratio)))))])
        disp(['Largest Reynolds number ' num2str(max(max(reynolds_check(~isnan(reynolds_check)))))])
        
        ratio_store_sec1(ratio_store_sec1>1.15) = NaN; %1.15
        
        %% Highlight 15 % Error Overshoot
        Q = bwperim(bwmorph(isnan( ratio_store_sec1' ),'fill')); % Using the contrast between the boundary to find the outline.
        Q(1,:) = 0;
        Q(end, :) = 0;
        Q(:,1) = 0;
        Q(:,end) = 0;
        limit_points = find(Q == 1);
        X = RHORHO(limit_points)./rho_0;
        Y = RR(limit_points);
        %         return
        ratio_store_sec2(ratio_store_sec2>1.15) = NaN; %1.15
        
        Q = bwperim(bwmorph(~isnan( ratio_store_sec2' ),'fill')); % Using the contrast between the boundary to find the outline.
        Q(1,:) = 0;
        Q(end, :) = 0;
        Q(:,1) = 0;
        Q(:,end) = 0;
        limit_points = find(Q == 1 & v_sec2 > 2);
        X_vc = RHORHO(limit_points)./rho_0;
        Y_vc = RR(limit_points);
        
        scatter([X; X_vc], [Y; Y_vc],[],[0.4940 0.1840 0.5560],'filled')
        
        %% Red-Dash Box for Commonly Used Particle Regions
        
        s= patch([1 1000 1000 1],[0 0 0.25 0.25],'red');
        s.FaceColor = 'none';
        s.LineStyle = '--';
        s.EdgeColor = 'red';
        s.LineWidth = 3;
        if jj == 2 && ii == 1
            interpn(r_normalised, rho_screen, Sc_size_sec2, 0.12, 19)
        end
        
        %% Labelling
        if jj ==1
            ylabel('r/\lambda [-]')
            yticks('auto')
            yticklabels('auto')
        else
            set(gca,'yticklabel',[])
        end
        if ii ==2
            xlabel('\rho_p / \rho_0 [-]')
            xticks('auto')
            xticklabels('auto')
            caxis([0 260])
            h123 = colorbar('Position', [0.6916+0.2134+0.01  0.1100+0.04  0.025  0.1100+0.2134*1.1], ...
                'YTick',[0:50:250],'FontSize',30);
            title(h123, 'S_r [mm]','FontSize',25);
            clabel(C,h, v_h,'FontSize',12,'Color','white','LabelSpacing',800)
        else
            set(gca,'xticklabel',[])
            caxis([0 120])
            h256 = colorbar('Position', [0.6916+0.2134+0.01  0.5200+0.05  0.025  0.1100+0.2134*1.1], ...
                'YTick',[0:20:120],'FontSize',30);
            title(h256, 'S_r [mm]','FontSize',25);
            clabel(C,h, v_h,'FontSize',12,'Color','white')
        end
        
        ylim([min(r_normalised) max(r_normalised)])
        xlim([min(rho_screen)./rho_0 max(rho_screen)./rho_0])
        grid on
        grid minor
        hold on
        set(gca,'layer','top')
        set(gca,'FontSize',30)
        
        if ii == 1
            title([num2str(frequency_sweep(jj)./1e03) ' kHz'],'FontSize',30)
        end
        if jj == 1
            text(0.01, 0.4, [num2str(pressure_sweep(ii)/1000) 'kPa'],'FontSize',30,'FontWeight','bold','Rotation',90)
        end
        
        tc_count = tc_count+1;
    end
end
%% Export Data
fig= gcf;
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[-0.125 -0.125 50 30],...
    'PaperSize',[50 30]);
print(fig,['results\main_figure\FIG4_APL_screen_multi_freq'],'-dpdf','-r0')
