clear all
close all
clc

% This files were tested on MATLAB2018a. 

% Toolbox Packages you may require:

% If testing analysis_file_generator;
% Symbolic Math Toolbox

% If reevaluating values calculated;
% Parallel Computing Toolbox (change parfor to for in the analysis file if you do not have it)

% Please contact Tatsuki FUSHIMI (tatsuki@levitation.engineer) for inquiry

load('air_properties_locked_data.mat','rho_0', 'v_sound', 'visc');

pressure_sweep = [5e03 25e03];      % Pressure amplitude in Pa
frequency_sweep = [20e03 40e03];    % Acoustic Frequency in Hz
r_min = 1/100;                      % min normalised radius
r_max = 1;                          % max nomalised radius
r_incr = 600;                       % number of increments in radius
rho_min = 1e-1;                     % min density in kg m3
rho_max = 1000;                     % max density in kg m3
rho_incr = 400;                     % number of increments in density
save_data = 1;                      % Export data to .mat file
%% Set Rendering Frequency
f_r = 10;                           % Rendering Frequency in Hz or FPS

%% Run Analysis
% Analysis files can be regenerated using:
% analysis_file_generator

% evaluation_function
%% Using data output, interpret data
% ----  Generate Fig. 1 to 3 ---- 
% Resultant graphs are in results/main_figure

graph_gen_FIG1_3

% ----  Generate Fig. 4 ---- 
% Requires tight_subplot.m from:
% Pekka Kumpulainen (2020). tight_subplot(Nh, Nw, gap, marg_h, marg_w) (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w), MATLAB Central File Exchange. Retrieved March 19, 2020.
% Download the .m file from above and put it in the same directory as this
% folder. 
% Resultant graphs are in results/main_figure

graph_gen_FIG4

% ----  Generate Fig. 5 ---- 
% Turn on analysis 

reanalyse_5 = 0; 
graph_gen_FIG5(reanalyse_5)

%% Generate supplementary figures
% Resultant graphs are in results/supp_figure
% ---- a) Calculation of the acoustic radiation force on a spherical particle ---- 

% supp_test_chen_calculation

% ---- b) Comparing rigid and elastic particle assumptions ---- 
% supp_test_rigid_elastic_check

% ---- c) Drag Coefficient ---- 
% supp_test_drag_calculation

% ---- e) Differences between vertical and horizontal rastering ---- 
% supp_test_compare_vertical_vs_horizontal

% ---- g) Extended parameter space to include \rho_p/\rho_0 \leq 1 ----
% supp_test_gen_density_sub1_fig1_3
% supp_test_gen_density_sub1_fig4

% ---- h) Effect of total sum of force exceeding the acoustic radiation force ---- 
% supp_recalculate = 0; % Set to 1 to recalculate
% supp_test_force_sepp(supp_recalculate)

% ---- i) Effect of Discretizing Screen Size ----
% supp_test_discretize

% ---- j) Analysing the cause of constant particle properties in Fig. 5 ----
% supp_test_lift_off

% ---- k) Case study when particle density is fixed at \rho_p = 19 kg m^{-3} ----
% supp_test_fixed_density

%% Additional Information
% Additional information about how to change the simulation file from Rigid
% assumption to non-rigid assumption (or generalised form) is available in

% additional_info.m

%% Identifying the file dependencies
% You can identify the file dependency of the matlab file using:
% [fList,pList] = matlab.codetools.requiredFilesAndProducts('XXXX.m');
% fList contains all function dependent