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

function [model] = sub_fig1_code()
%% This function generates data points for zero ARF points for p_a = 5 kPa, and f_f = 40 kHz. 
p_a = 5000;
acoustic_f = 40e03;

%% Analyse Results
load(['results/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa.mat']) %#ok<LOAD>
lambda = (v_sound/acoustic_f);
r = linspace(lambda/100, lambda,700);
rhorho_zero_arf = [];
rr_zero_arf = [];
epsilon = 1e-10;
for ii = 1:length(rho_screen)
    rho_p = rho_screen(ii);
    F_arf = abs(alpha_c(m_sum_calculator(r, acoustic_f, v_sound, rho_0, rho_p, epsilon),p_a,rho_0, acoustic_f));
    [pks,locs] = findpeaks(-F_arf); %Find minimas in the ARF along r/lambda
    if ~isempty(locs)
        locs((-pks)>0.01e-03) = [];
        rr_zero_arf = [rr_zero_arf r(locs)./lambda]; %#ok<AGROW>
        rhorho_zero_arf = [rhorho_zero_arf rho_p.*ones(1,length(locs))./rho_0]; %#ok<AGROW>
    end
end
%% output function
%% Chopped the lines into sectors such that a nice line can be drawn.
model = struct;
%% Section 0.85 < r/lambda 
idx = find(rr_zero_arf > 0.85);
model.sect1_rr = rr_zero_arf(idx);
model.sect1_rho = rhorho_zero_arf(idx);
%% Section 2 0.7 <= r/lambda < 0.85
idx = find(rr_zero_arf >= 0.7 & rr_zero_arf < 0.85);
model.sect2_rr = rr_zero_arf(idx);
model.sect2_rho = rhorho_zero_arf(idx);
%% Section 3 0.6 < r/lambda <= 0.7
idx = find(rr_zero_arf > 0.6 & rr_zero_arf <= 0.702);
model.sect3_rr = rr_zero_arf(idx);
model.sect3_rho = rhorho_zero_arf(idx);
%% Section 4 0.38 < r/lambda < 0.5
idx = find(rr_zero_arf > 0.38 & rr_zero_arf < 0.5);
model.sect4_rr = rr_zero_arf(idx);
model.sect4_rho = rhorho_zero_arf(idx);
%% Section 5 r/lambda < 0.35
idx = find(rr_zero_arf < 0.35);
model.sect5_rr = rr_zero_arf(idx);
model.sect5_rho = rhorho_zero_arf(idx);
end

