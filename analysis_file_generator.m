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

clear all
close all
clc

%% Test Condition: MATLAB R2018a, WIN 10
% Using Symbolic Math Toolbox version 8.1.

syms rho_p screen_freq screen_pixel t acoustic_f rho_0 visc v_sound F2 m p_a a_z vel_z r

k = 2*pi*acoustic_f./v_sound;
lambda = (v_sound/acoustic_f);
volume = (4/3)*pi*((r)^3);
mass =  volume*rho_p;
virtual_mass = 0.5*volume*rho_0; %Half the volume of the body, resulting from the surrounding fluids (i.e. air)
%% Acoustic Radiation Force
% Using Chen & Apfel (Chen, X. & Apfel, R. E. Radiation force on a spherical object in an axisymmetric wave field and its application to the calibration of high-frequency transducers. J. Acoust. Soc. Am. 99, 713 (1996).) 
% F2 is the 'infinite sum' in Chen & Apfel formulation.
alpha_c = (p_a^2).*F2.*pi ./ (rho_0.*((2*pi*acoustic_f).^2)); 

%% Gravity Force + Buyoncy Force (Balancing gravitational force and buoyancy force) 
% Worst case scenario can either happen when going up (buoynacy rho_0 > rho_p) or down (gravity (rho_p > rho_0)). 
% Thus, obtained the negative of the total amplitude of the two forces combined. 
% Using formulation by: McKee, K. & Czarnecki, A. Acceleration due to buoyancy and mass renormalization. Am. J. Phys. 87, 165–170 (2019).
g = 9.81; 
F_gravity_buoyancy = -abs((-mass + volume*rho_0)*g); 
%% Drag Force Calculation
% Coeffcient of drag from: Morrison, F. A. Data correlation for drag coefficient for sphere. Department of Chemical Engineering, Michigan Technological University, Houghton, MI 49931, 1–2 (2013).
% Use test_drag_calculation.m for comparison against Flemmer & Banks (Flemmer, R. L. C. & Banks, C. L. On the drag coefficient of a sphere. Powder Technol. 48, 217–221 (1986).)
Re = (vel_z .* rho_0 .* (2.*r) ./(visc));
C_d = 24./(Re) + (2.6.*Re./5.0)./(1+(Re./5.0).^1.52) +...
    (0.411.*((Re)./(2.63e05)).^-7.94)./(1+(Re./(2.63e05)).^-8) + ...
    (0.25.*(Re./1e06))/(1+(Re./1e06));
F_drag = (-1 .* C_d .* (pi/4) .* (2.*r).^2 .* 0.5 .* rho_0 .* (vel_z.^2));

Re_abs = abs(vel_z .* rho_0 .* (2.*r) ./(visc));
C_d_abs = 24./(Re_abs) + (2.6.*Re_abs./5.0)./(1+(Re_abs./5.0).^1.52) +...
    (0.411.*((Re_abs)./(2.63e05)).^-7.94)./(1+(Re_abs./(2.63e05)).^-8) + ...
    (0.25.*(Re_abs./1e06))/(1+(Re_abs./1e06));
F_drag_abs = (-1 .* C_d_abs .* (pi/4) .* (2.*r).^2 .* 0.5 .* rho_0 .* (vel_z.^2).*sign(vel_z));


%% Generate Function Files
matlabFunction(alpha_c,'File','alpha_c','Vars',[F2,p_a,rho_0,acoustic_f]);
matlabFunction(F_drag,'File','F_drag','Vars',[r,rho_0,vel_z,visc]);
matlabFunction(F_drag_abs,'File','F_drag_abs','Vars',[r,rho_0,vel_z,visc]);
matlabFunction(F_gravity_buoyancy,'File','F_gravity_buoyancy','Vars',[r,rho_p, rho_0]);
matlabFunction(Re,'File','Reynolds','Vars',[r,rho_0,vel_z,visc]);
matlabFunction(virtual_mass,'File','virtual_mass','Vars',[r,rho_0]);
matlabFunction(mass,'File','mass','Vars',[r,rho_p]);

clear all 
syms Re

C_d = 24./(Re) + (2.6.*Re./5.0)./(1+(Re./5.0).^1.52) +...
    (0.411.*((Re)./(2.63e05)).^-7.94)./(1+(Re./(2.63e05)).^-8) + ...
    (0.25.*(Re./1e06))/(1+(Re./1e06));

matlabFunction(C_d,'File','C_d','Vars',[Re]);