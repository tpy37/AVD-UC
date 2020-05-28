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
%% Evaluates the discretized version of the screen size
load('air_properties_locked_data.mat','rho_0', 'v_sound', 'visc');

pressure_sweep = [5e03 25e03];      % Pressure amplitude in Pa
frequency_sweep = [20e03 40e03];    % Acoustic Frequency in Hz
%% Set Rendering Frequency
f_r = 10;

tc_count = 1;
figure('units','normalized','outerposition',[0 0 1 1])
ha = tight_subplot(2,2,[.04 .04],[.15 .08],[.1 .12]);

v_h = [0:40:260];
cyc_len = 300;

for ii = 1:length(pressure_sweep)
    p_a = pressure_sweep(ii);
    for jj = 1:length(frequency_sweep)
        acoustic_f = frequency_sweep(jj);
        load(['results/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa.mat'])
        
        if tc_count == 1
            [RR, RHORHO]= ndgrid(r_normalised, rho_screen);
        end
        
        wavelength = v_sound/acoustic_f;
        
        %% Sec 1
        Sc_size_sec1 = 2.*sqrt((max_velocity_store'.*(RR.*wavelength))./(pi.*f_r));
        omega_sec1 = sqrt((f_r.*pi.*max_velocity_store')./(RR.*(wavelength)));
        Sc_size_sec1(omega_sec1>omega_n_store') = NaN;
        Sc_size_sec1(RHORHO<1) = NaN;
        
        
        %% Sec 2
        v_sec2 = ((((max_velocity_store'.*omega_n_store').^2).*RR.*wavelength)./(f_r.*pi)).^(1/3);
        v_sec2(v_sec2 > max_velocity_store') = NaN;
        omega_ac = (1./(v_sec2).*max_velocity_store'.*omega_n_store');
        Sc_size_sec2 = 2.*sqrt((v_sec2.*(RR.*wavelength))./(pi.*f_r));
        Sc_size_sec2(RHORHO<1) = NaN;
        
        %% N_par must be integer value. Recalculate screen size after making it into interger.
        N_par_sec1 = floor(Sc_size_sec1./(2.*RR.*wavelength));
        N_par_sec2 = floor(Sc_size_sec2./(2.*RR.*wavelength));
        
        S_sec1 = 2.*RR.*wavelength.*N_par_sec1;
        S_sec2 = 2.*RR.*wavelength.*N_par_sec2;
        %% Calculate the difference between the maximum screen size in continuous and discrete
        rat_s1 = max(max(Sc_size_sec1))./max(max(S_sec1));
        rat_s2 = max(max(Sc_size_sec2))./max(max(S_sec2));
        max(max(rat_s1(rat_s1 ~= Inf)))
        max(max(rat_s2(rat_s2 ~= Inf)))
        
        %% Plot
        axes(ha(tc_count))
        hold on
        [C, h] = contourf(rho_screen./rho_0, r_normalised, S_sec1.*1000,[0:10:260]);
        [Cv, hv] = contourf(rho_screen./rho_0, r_normalised, S_sec2.*1000,[0:10:260]);
        hold on
        shading interp
        set(gca, 'XScale', 'log')
        
        if jj == 2 && ii == 1
            interpn(r_normalised, rho_screen, Sc_size_sec2, 0.12, 19);
        end
        
        if jj ==1
            ylabel('r/\lambda [-]')
            yticks('auto')
            yticklabels('auto')
        else
            set(gca,'yticklabel',[])
        end
        
        if or(tc_count==1, tc_count == 2)
            caxis([0 60])
            clabel(Cv,hv, v_h,'FontSize',12,'Color','white','LabelSpacing',300)
            v_h = [0:20:120];
            set(gca,'xticklabel',[])
        end
        
        if or(tc_count == 3, tc_count ==4)
            caxis([0 180])
            xlabel('\rho_p / \rho_0 [-]')
            xticks('auto')
            xticklabels('auto')
            v_h = [0:40:180];
            clabel(Cv,hv, v_h,'FontSize',12,'Color','white','LabelSpacing',9000)
        end
        
        if tc_count == 2
            h123 = colorbar('Position', [0.6916+0.2134+0.01  0.5200+0.05  0.025  0.1100+0.2134*1.1], ...
                'YTick',[0:10:60],'FontSize',30);
            title(h123, 'S_r [mm]','FontSize',25);
        end
        
        if tc_count == 4
            h256 = colorbar('Position', [0.6916+0.2134+0.01  0.1100+0.04  0.025  0.1100+0.2134*1.1], ...
                'YTick',[0:20:180],'FontSize',30);
            title(h256, 'S_r [mm]','FontSize',25);
        end
        
        ylim([min(r_normalised) max(r_normalised)])
        xlim([1 max(rho_screen)./rho_0])
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

saveas(gcf,'results\supp_figure\discretize.png')
