%% run_postprocess.m
% Quick-reference script for all post-processing functions.
% Run this after completing simulations.
%
% Each section can be run independently (select and Ctrl+Enter).

addpath(genpath('./Models'))
addpath(genpath('./Shapes'))
addpath(genpath('./PostProcessing'))

opts = pub_style();

%% ======================================================
%  1. Single-scenario full post-processing
% ======================================================
% Generates Fig 01-06 in ./Figures/
postprocess_seal('./Results_Seal/Standard_21MPa_80C/');


%% ======================================================
%  2. Custom evolution snapshots (specify time steps)
% ======================================================
result_folder = './Results_Seal/Standard_21MPa_80C/';
steps = [10 50 100 200];  % Which saved steps to show
plot_evolution_snapshots(result_folder, steps, opts, './Figures');


%% ======================================================
%  3. Line profiles across the seal
% ======================================================
data = load('./Results_Seal/Standard_21MPa_80C/end.mat');
plot_line_profiles(data, opts, './Figures');


%% ======================================================
%  4. Multi-scenario comparison
% ======================================================
folders = {'./Results_Seal/Standard_21MPa_80C/', ...
           './Results_Seal/HighTemp_21MPa_120C/', ...
           './Results_Seal/HighPress_35MPa_80C/'};
labels = {'Standard', 'High-Temp', 'High-Press'};
compare_scenarios(folders, labels, opts);


%% ======================================================
%  5. LaTeX results table
% ======================================================
labels_tex = {'Standard (21 MPa, 80$^\circ$C)', ...
              'High-temp (21 MPa, 120$^\circ$C)', ...
              'High-press (35 MPa, 80$^\circ$C)'};
generate_results_table(folders, labels_tex, './Figures/table_results.tex');


%% ======================================================
%  6. Animation (MP4 video)
% ======================================================
% Phase-field evolution animation
generate_animation('./Results_Seal/Standard_21MPa_80C/', 'phi', ...
    'FrameRate', 12, 'StepInterval', 2);

% Fluid diffusion animation
generate_animation('./Results_Seal/Standard_21MPa_80C/', 'CL', ...
    'FrameRate', 12, 'StepInterval', 2);

% Aging animation
generate_animation('./Results_Seal/Standard_21MPa_80C/', 'alpha_a', ...
    'FrameRate', 12, 'StepInterval', 5);


%% ======================================================
%  7. Custom single-panel figures (for slides / posters)
% ======================================================
data = load('./Results_Seal/Standard_21MPa_80C/end.mat');

% --- Large single damage field ---
fig = figure; set(fig, 'Position', opts.fig_square, 'Color', 'w');
plot_nodal_custom(data.physics, data.mesh, "phi", 0, "Internal", opts.m2mm);
colormap(flipud(hot(256))); caxis([0 1]);
cb = pub_colorbar('$\phi$ [-]', opts);
xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, opts.axes_props{:}); axis image;
export_fig(fig, 'single_phi', opts, 'both');

% --- Large deformed shape ---
fig = figure; set(fig, 'Position', opts.fig_square, 'Color', 'w');
plot_nodal_custom(data.physics, data.mesh, "phi", 50, "Internal", opts.m2mm);
colormap(flipud(hot(256))); caxis([0 1]);
cb = pub_colorbar('$\phi$ [-]', opts);
xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, opts.axes_props{:}); axis image;
title('Deformed ($\times 50$)', 'Interpreter', 'latex', 'FontSize', 14);
export_fig(fig, 'single_deformed', opts, 'both');


%% ======================================================
%  8. Export settings overview
% ======================================================
%  Format options for export_fig:
%    'raster'  - PNG at 600 DPI (default, for contour plots)
%    'vector'  - EPS (for line plots, editable in Illustrator)
%    'pdf'     - PDF with best-fit
%    'both'    - PNG + EPS simultaneously
%    'tiff'    - TIFF at 600 DPI (some journals require this)
%
%  To change DPI:
%    opts.dpi = 300;   % Lower res for drafts
%    opts.dpi = 1200;  % Ultra-high for print
%
%  To change font:
%    opts.fontname = 'Arial';        % Sans-serif alternative
%    opts.fontname = 'Helvetica';    % Another sans option


function plot_nodal_custom(physics, mesh, dofName, dispscale, plotloc, coord_scale)
    [dxTypes, dxSteps] = physics.dofSpace.getDofType({"dx";"dy";dofName});
    for g = 1:length(mesh.Elementgroups)
        if mesh.Elementgroups{g}.name == plotloc && mesh.Elementgroups{g}.type == "Q9"
            for el = 1:size(mesh.Elementgroups{g}.Elems, 1)
                elnodes = mesh.Elementgroups{g}.Elems(el,:);
                order = [1 3 9 7];
                zdofs = physics.dofSpace.getDofIndices(dxTypes(3), elnodes);
                if dispscale >= 0
                    xdofs = physics.dofSpace.getDofIndices(dxTypes(1), elnodes);
                    ydofs = physics.dofSpace.getDofIndices(dxTypes(2), elnodes);
                    X(el,:) = (mesh.Nodes(elnodes(order),1) + ...
                        dispscale*physics.StateVec{dxSteps(1)}(xdofs(order))) * coord_scale;
                    Y(el,:) = (mesh.Nodes(elnodes(order),2) + ...
                        dispscale*physics.StateVec{dxSteps(2)}(ydofs(order))) * coord_scale;
                else
                    X(el,:) = mesh.Nodes(elnodes(order),1) * coord_scale;
                    Y(el,:) = mesh.Nodes(elnodes(order),2) * coord_scale;
                end
                Z(el,:) = physics.StateVec{dxSteps(3)}(zdofs(order));
            end
            patch(X', Y', Z', Z', 'EdgeColor', 'None', 'FaceColor', 'interp');
            hold on;
        end
    end
end
