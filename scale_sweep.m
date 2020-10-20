clearvars;
close all hidden;
clc;

format shortEng;

%% plotting settings:
addpath("\\home.ansatt.ntnu.no\geraldod\Documents\MATLAB\scripts");
plot_definitions;
fig_dim = figure_dimensions("custom", "landscape", 16.0);

clr = linspecer(6, 'qualitative');
plot_prop1 = plot_properties("line", 1, 1, line_style, 2.0, clr);
plot_prop2 = plot_properties("line", 2, 2, line_style, 2.0, clr);

%% scaling:
ref_DT = NREL_5MW();

gamma_P = scaling_factor.generate_gamma(4, 100);
gamma_P = fliplr(gamma_P);
P_scale = gamma_P*ref_DT.P_rated;

[gamma, res, SH, f_n, mode_shape, k_mesh, gamma_asp, obj_array] = ...
    ref_DT.scaled_sweep(P_scale,'aspect_set', ["integrity", ...
                                               "integrity_refine", ...
                                               "resonances", ...
                                               "resonances_refine"]);
                                       
save('data_scale_sweep_LP', 'gamma', 'res', 'SH', 'f_n', ...
    'mode_shape', 'k_mesh', 'gamma_asp', 'obj_array', 'P_scale');