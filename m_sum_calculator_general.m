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

function [F2] = m_sum_calculator_general(r, acoustic_f, v_sound, rho_0, rho_p, c_c, c_s, epsilon)
%% Calculates infintie sum in Chen & Apfel till it reaches convergence rate of epsilon. 
%% For non-rigid assumption. 
%% Inputs
% r - radius (in m)
% acoustic_f - acoustic frequency (in Hz)
% v_sound - speed of sound in fluid (in m s^{-1})
% rho_0 - density of fluid (in kg m^{-3})
% rho_p - density of particle (in kg m^{-3})
% c_c - comparessional wave speed in particle (in m s^{-1})
% c_s - shear wave speed in particle (in m s^{-1})
% epsilon - termination point for infinite sum


%% Outputs
% F2 - The converged infinite sum in Chen & Apfel.
%% Initialize parameters
m_sum = 0;
convergence_rate = 1;
m = 0;
while convergence_rate > epsilon
    % Cacculate D_m and D_{m+1}
    Dm = calculate_D_general(2*pi*acoustic_f./v_sound,r, rho_0, rho_p, 2*pi*acoustic_f./c_c, 2*pi*acoustic_f./c_s, m);
    Dm1 = calculate_D_general(2*pi*acoustic_f./v_sound,r, rho_0, rho_p, 2*pi*acoustic_f./c_c, 2*pi*acoustic_f./c_s, m+1);
    
    % Calculate alpha_m and beta_m
    Bm = 1i.*Dm./(1-1i.*Dm);
    alpha_m = real(Bm);
    beta_m = imag(Bm);
    
    % Calculate alpha_{m+1) and beta_{m+1)
    Bm1 = 1i.*Dm1./(1-1i.*Dm1);
    alpha_m1 = real(Bm1);
    beta_m1 = imag(Bm1);
    
    % Store previous infinite sum m_sum
    prev_m_sum = m_sum;
    % Update m_sum
    m_sum = m_sum + (m+1).*((-1).^m).*(-beta_m+beta_m1+2.*alpha_m.*beta_m1-2.*beta_m.*alpha_m1);
    
    % Calculate convergence rate, if max of vector is abov epsilon,
    % continue. 
    convergence_rate = max(abs(prev_m_sum-m_sum)./(abs(m_sum)));
    m = m+1;
    
    % To visualize, uncomment.
%     figure(3)
%     hold on
%     scatter(m, m_sum)
end
%Output the infinite sum
F2 = m_sum;
end