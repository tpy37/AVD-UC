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

function [omega_n_store, max_velocity_store] = main_function_horizontal(p_a,acoustic_f, rho_0, v_sound, visc, r_min, r_max, r_incr, rho_min, rho_max, rho_incr, save_data)
%% Main Function for paper, "Theoretical Consideration towards the Designs of Acoustohporetic Volumetric Display"
%% This function evaluates the maximum achievable velocity, generalised frequency, and generalised amplitude based on the assumptions
%% as stated in the original manuscripts.
%% Takes input, p_a (Pressure Amplitude in Pa), and acoustic_f (Acoustic Frequency in Hz),
%% rho_0 (density in kg m^{-3}), v_sound (speed of sound in ms ^{-1}), visc (viscosity in kg (m s)^{-1})
%% Returns max_velocity_store (maximum achievable velocity) and 'omega_n_store' (maximum achievable generalised frequency)
%% Also exports all data to (.mat) file if 'save_data = 1'.
%% M_SUM TOLERANCE
m_sum_tol = 1e-10;

%% Newton-Raphson Settings
iteration_max = 5000;
lit_h =0.05;
rel_v_d = 1e-10;

%% Radius Search Range
r_normalised = linspace(r_min, r_max, r_incr);

%% Density Search Rnage
rho_screen = logspace(log10(rho_min), log10(rho_max), rho_incr);

%% Acoustic Parameters
wavelength = v_sound/acoustic_f; % wavelength in m

%% Non-normalised radius
r = r_normalised.*wavelength;

%% Data Export
max_velocity_store = zeros(rho_incr, r_incr); % Stores maximum achievable velocity
omega_n_store = zeros(rho_incr, r_incr); % Stores maximum achievable generalised frequency

arf_store = zeros(rho_incr, r_incr); % Stores acoustic radiation force at the point
drag_store = zeros(rho_incr, r_incr); % Stores the drag force at the point
reynolds_check = zeros(rho_incr, r_incr); % Stores the value of Reynolds number at the point
rem_check = zeros(rho_incr, r_incr); % Stores the functional value at the optimised point

for rh = 1:length(rho_screen) % Sweep in rho
    rho_p = rho_screen(rh);
    parfor ii = 1:length(r) % Sweep in radius
        
        %% Evaluate M Sum for ARF
        [F2] = m_sum_calculator(r(ii), acoustic_f, v_sound, rho_0, rho_p, m_sum_tol);
        
        %% Record ARF data
        arf_store(rh, ii) = abs(alpha_c(F2,p_a,rho_0,acoustic_f));  %#ok<PFOUS>
        
        %% Setup functions for Newton Raphson Optimization
        function_t = @(v_cand) abs(alpha_c(F2,p_a,rho_0,acoustic_f))+Drag_Force_Zero_Escape(1, r(ii),rho_0,v_cand,visc);
        derivative_t = @(v_cand) (function_t(v_cand+lit_h) - function_t(v_cand))/(lit_h);
        
        %% Reset counter of the optisation to 1, and set initial guess of velocity as 1 m/s
        counter = 1;
        orig_v = 15;
        v_cand = orig_v;
        
        %% Loop until convers or reaches maxima.
        while 1
            prev_v = v_cand;
            %% Newton Raphson Function
            v_cand = v_cand - function_t(v_cand)/derivative_t(v_cand);
            
            %% "NaN" out, if invalid solutions (i.e. velocity < 0 )
            if v_cand < 0
                v_cand = NaN;
                break
            end
            
            %% If convergence rate is below tolerance, exit
            if abs(prev_v-v_cand)/(v_cand) < rel_v_d
                break
            elseif counter> iteration_max
                counter = 1;
                orig_v = orig_v + 10;
                v_cand = orig_v;
            end
            %% Counter adds
            counter = counter + 1;
        end
        
        %% Record Data
        rem_check(rh, ii) = function_t(v_cand);  %#ok<PFOUS>
        drag_store(rh, ii) = Drag_Force_Zero_Escape(1, r(ii),rho_0,v_cand,visc);  %#ok<PFOUS>
        max_velocity_store(rh, ii) = v_cand;
        reynolds_check(rh, ii) = Reynolds(r(ii),rho_0,v_cand,visc);
        
        %% Evaluate generalised frequency omega_g using maximum achievable velocity
        max_alpha = abs(alpha_c(F2,p_a,rho_0,acoustic_f));
        Mass_total = mass(r(ii),rho_p)+virtual_mass(r(ii),rho_0); % inertia needs to consider both particle and added mass
        omega_n_max = max_alpha./(Mass_total*v_cand);
        %% Store maximum achievable vleocity
        omega_n_store(rh, ii) = omega_n_max;
    end
end

%% Evaluate generalised amplitude A using maximum achievable velocity and generalised frequency
A_pix = max_velocity_store./omega_n_store;  %#ok<NASGU>
%% If maximum value of reynolds number exceeds the valid range of the equation, return error.
if max(max(reynolds_check)) > 1e06
    error('Max Reynolds ERROR')
elseif max(max(max_velocity_store)) >= v_sound
    warning('Max Velocity Exceeds Speed of Sound. Check Solution')
else
    %% Export data into study 1
    if save_data==1
        [status, ~, ~] = mkdir('results');
        if status
            save(['results/study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa_hoz.mat'], 'A_pix', 'omega_n_store',...
                'max_velocity_store','drag_store','arf_store',...
                'reynolds_check','rho_screen','r','rem_check','r_min','r_max', ...
                'r_incr', 'rho_min', 'rho_max','rho_incr','r_normalised','rho_0',...
                'v_sound', 'visc');
            disp(['study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa'])
            disp('----- Calculation Completed and Exported without Error -----')
        else
            disp('Saving Error')
        end
    else
        disp(['study_ff_' num2str(acoustic_f) '_Hz_pp_' num2str(p_a) '_Pa'])
        disp('----- Calculation Completed without Error (No .mat File generated, save_data = 1 for generation) -----')
    end
end

end

