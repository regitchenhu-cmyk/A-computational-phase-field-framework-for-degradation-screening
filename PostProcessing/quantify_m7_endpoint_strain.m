function summary = quantify_m7_endpoint_strain(outputDir)
%QUANTIFY_M7_ENDPOINT_STRAIN Recover endpoint strains without rerunning cases.
%
% This post-processor addresses Major Comment 7 by reading the saved
% end.mat files and evaluating the displacement gradient at every interior
% Gauss point.  It reports a robust whole-domain P99 together with the raw
% point maximum and its coordinates, so a constrained-corner hotspot is not
% confused with a representative strain level.
%
% The reported infinitesimal von Mises-equivalent strain is
%
%   eps_eq = sqrt((2/3) * dev(eps):dev(eps)),
%
% where eps = sym(grad(u)) and eps_zz = 0 (plane strain).  Isotropic
% swelling does not change this deviatoric equivalent strain.  The script
% also computes Green-Lagrange strains from F = I + grad(u) as a kinematic
% comparison; this is a post-processing diagnostic, not a rerun or a
% finite-strain validation of the coupled model.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(repoRoot, ...
        'Results_Paper_ScientificallyCorrected_20260823', ...
        'M7_Strain_Quantification_20260828');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'Models')));
addpath(genpath(fullfile(repoRoot, 'Shapes')));

studyRoots = {
    'MainCampaign', fullfile(repoRoot, ...
        'Results_Parametric_ScientificallyCorrected_20260823');
    'PaperChecks', fullfile(repoRoot, ...
        'Results_Paper_ScientificallyCorrected_20260823')
};

caseRows = struct([]);
referencePoints = table();
for r = 1:size(studyRoots, 1)
    studySet = string(studyRoots{r, 1});
    studyRoot = studyRoots{r, 2};
    files = dir(fullfile(studyRoot, '**', 'end.mat'));
    for k = 1:numel(files)
        endFile = fullfile(files(k).folder, files(k).name);
        relCase = erase(string(files(k).folder), string(studyRoot) + filesep);
        if relCase == string(files(k).folder)
            relCase = string(files(k).folder);
        end
        fprintf('M7 strain recovery: %s / %s\n', studySet, relCase);
        [row, pointTable] = analyze_endpoint(endFile, studySet, relCase);
        if isempty(caseRows)
            caseRows = row;
        else
            caseRows(end + 1, 1) = row; %#ok<AGROW>
        end
        if studySet == "MainCampaign" && relCase == "Grid_G28x40"
            referencePoints = pointTable;
        end
    end
end

if isempty(caseRows)
    error('M7Strain:NoEndpoints', 'No end.mat files were found.');
end

summary = struct2table(caseRows);
summary = sortrows(summary, {'StudySet', 'Case'});
writetable(summary, fullfile(outputDir, 'm7_strain_summary.csv'));

meshMask = summary.StudySet == "MainCampaign" & ...
    startsWith(summary.Case, "Grid_G");
meshSummary = summary(meshMask, :);
writetable(meshSummary, fullfile(outputDir, 'm7_mesh_strain_summary.csv'));

if ~isempty(referencePoints)
    writetable(referencePoints, ...
        fullfile(outputDir, 'm7_reference_Grid_G28x40_ip_strain.csv'));
end

save(fullfile(outputDir, 'm7_strain_quantification.mat'), ...
    'summary', 'meshSummary', 'referencePoints');

write_readme(outputDir, summary, meshSummary);
disp(summary(:, {'StudySet', 'Case', 'FinalTime_days', ...
    'EqvP99_pct', 'EqvMax_pct', 'EqvMaxX_mm', 'EqvMaxY_mm', ...
    'FracEqvGt30_pct', 'GreenVsSmallP99_pct'}));
end

function [row, pointTable] = analyze_endpoint(endFile, studySet, caseName)
s = load(endFile, 'mesh', 'physics');
mesh = s.mesh;
physics = s.physics;
model = find_model(physics, "SealPhaseFieldDamage");
if isempty(model)
    error('M7Strain:MissingModel', ...
        'SealPhaseFieldDamage was not found in %s.', endFile);
end

group = model.myGroupIndex;
elems = mesh.Elementgroups{group}.Elems;
nElem = size(elems, 1);
nIp = mesh.Elementgroups{group}.ShapeFunc.ipcount;
nPoint = nElem * nIp;

[uTypes, uSteps] = physics.dofSpace.getDofType({'dx', 'dy'});
if any(uTypes == 0)
    error('M7Strain:MissingDofs', 'Required displacement dofs are absent.');
end
if uSteps(1) ~= uSteps(2)
    error('M7Strain:StepMismatch', 'dx and dy do not share a solution step.');
end

x = zeros(nPoint, 1);
y = zeros(nPoint, 1);
element = zeros(nPoint, 1);
ipLocal = zeros(nPoint, 1);
eqvSmall = zeros(nPoint, 1);
eqvGreen = zeros(nPoint, 1);
maxAbsPrincipal = zeros(nPoint, 1);
greenVsSmall = zeros(nPoint, 1);
detF = zeros(nPoint, 1);
rotationNorm = zeros(nPoint, 1);

cursor = 0;
for e = 1:nElem
    nodes = mesh.getNodes(group, e);
    [~, G, w] = mesh.getVals(group, e);
    coords = mesh.getIPCoords(group, e);
    uxDofs = physics.dofSpace.getDofIndices(uTypes(1), nodes);
    uyDofs = physics.dofSpace.getDofIndices(uTypes(2), nodes);
    ux = physics.StateVec{uSteps(1)}(uxDofs);
    uy = physics.StateVec{uSteps(1)}(uyDofs);

    for ip = 1:numel(w)
        cursor = cursor + 1;
        dNdx = squeeze(G(ip, :, 1));
        dNdy = squeeze(G(ip, :, 2));
        dNdx = dNdx(:);
        dNdy = dNdy(:);

        gradU = [ux(:)' * dNdx, ux(:)' * dNdy; ...
                 uy(:)' * dNdx, uy(:)' * dNdy];
        eps2 = 0.5 * (gradU + gradU');
        eps3 = zeros(3);
        eps3(1:2, 1:2) = eps2;

        F = eye(3);
        F(1:2, 1:2) = eye(2) + gradU;
        Egreen = 0.5 * (F' * F - eye(3));

        eqvSmall(cursor) = equivalent_strain(eps3);
        eqvGreen(cursor) = equivalent_strain(Egreen);
        principal = eig(eps3);
        maxAbsPrincipal(cursor) = max(abs(principal));
        greenVsSmall(cursor) = 100 * norm(Egreen - eps3, 'fro') / ...
            max(norm(Egreen, 'fro'), 1e-14);
        detF(cursor) = det(F);
        rotationNorm(cursor) = norm(0.5 * (gradU - gradU'), 'fro');

        x(cursor) = coords(1, ip);
        y(cursor) = coords(2, ip);
        element(cursor) = e;
        ipLocal(cursor) = ip;
    end
end

if cursor ~= nPoint
    error('M7Strain:PointCount', 'Recovered point count is inconsistent.');
end

[eqvMax, iMax] = max(eqvSmall);
positiveJ = detF > 0;
if any(positiveJ)
    positiveIndices = find(positiveJ);
    [eqvMaxPositiveJ, localPositiveMax] = max(eqvSmall(positiveJ));
    iMaxPositiveJ = positiveIndices(localPositiveMax);
else
    eqvMaxPositiveJ = NaN;
    iMaxPositiveJ = NaN;
end
cornerMask = x <= 0.20e-3 & y <= 0.20e-3;
if ~any(cornerMask)
    error('M7Strain:CornerEmpty', ...
        'No integration points were found in the 0.20 mm corner window.');
end

row = struct();
row.StudySet = studySet;
row.Case = caseName;
row.SourceEndMat = portable_source_label(studySet, caseName);
row.FinalTime_days = physics.time / 86400;
row.Elements = nElem;
row.IntegrationPoints = nPoint;
row.EqvP95_pct = 100 * percentile(eqvSmall, 95);
row.EqvP99_pct = 100 * percentile(eqvSmall, 99);
row.EqvP999_pct = 100 * percentile(eqvSmall, 99.9);
row.EqvMax_pct = 100 * eqvMax;
row.EqvMaxX_mm = 1e3 * x(iMax);
row.EqvMaxY_mm = 1e3 * y(iMax);
row.EqvMaxElement = element(iMax);
row.EqvMaxIp = ipLocal(iMax);
row.EqvMaxNearestCorner = nearest_corner(mesh, x(iMax), y(iMax));
row.DetFAtEqvMax = detF(iMax);
row.EqvMaxPositiveJ_pct = 100 * eqvMaxPositiveJ;
if isfinite(iMaxPositiveJ)
    row.EqvMaxPositiveJX_mm = 1e3 * x(iMaxPositiveJ);
    row.EqvMaxPositiveJY_mm = 1e3 * y(iMaxPositiveJ);
else
    row.EqvMaxPositiveJX_mm = NaN;
    row.EqvMaxPositiveJY_mm = NaN;
end
row.CornerPointCount = nnz(cornerMask);
row.CornerEqvP99_pct = 100 * percentile(eqvSmall(cornerMask), 99);
row.CornerEqvMax_pct = 100 * max(eqvSmall(cornerMask));
row.MaxAbsPrincipalP99_pct = 100 * percentile(maxAbsPrincipal, 99);
row.MaxAbsPrincipalMax_pct = 100 * max(maxAbsPrincipal);
row.GreenEqvP99_pct = 100 * percentile(eqvGreen, 99);
row.GreenEqvMax_pct = 100 * max(eqvGreen);
row.GreenVsSmallP99_pct = percentile(greenVsSmall, 99);
row.GreenVsSmallAtEqvMax_pct = greenVsSmall(iMax);
row.DetFMin = min(detF);
row.DetFP01 = percentile(detF, 1);
row.DetFP99 = percentile(detF, 99);
row.DetFMax = max(detF);
row.NonpositiveDetFCount = nnz(~positiveJ);
row.NonpositiveDetF_pct = 100 * mean(~positiveJ);
row.RotationNormP99 = percentile(rotationNorm, 99);
row.FracEqvGt10_pct = 100 * mean(eqvSmall > 0.10);
row.FracEqvGt20_pct = 100 * mean(eqvSmall > 0.20);
row.FracEqvGt30_pct = 100 * mean(eqvSmall > 0.30);

pointTable = table(element, ipLocal, 1e3*x, 1e3*y, ...
    100*eqvSmall, 100*eqvGreen, 100*maxAbsPrincipal, ...
    greenVsSmall, detF, rotationNorm, cornerMask, ~positiveJ, ...
    'VariableNames', {'Element', 'IntegrationPoint', 'X_mm', 'Y_mm', ...
    'EqvSmall_pct', 'EqvGreen_pct', 'MaxAbsPrincipal_pct', ...
    'GreenVsSmall_pct', 'DetF', 'RotationNorm', 'InCorner0p2mm', ...
    'KinematicInversion'});
end

function value = equivalent_strain(epsTensor)
dev = epsTensor - trace(epsTensor) / 3 * eye(3);
value = sqrt((2/3) * sum(dev(:).^2));
end

function mdl = find_model(physics, modelName)
mdl = [];
for k = 1:numel(physics.models)
    candidate = physics.models{k};
    if isprop(candidate, 'myName') && string(candidate.myName) == modelName
        mdl = candidate;
        return;
    end
end
end

function q = percentile(values, pct)
values = values(isfinite(values));
if isempty(values)
    q = NaN;
    return;
end
q = prctile(values, pct);
end

function label = nearest_corner(mesh, x, y)
xmin = min(mesh.Nodes(:, 1));
xmax = max(mesh.Nodes(:, 1));
ymin = min(mesh.Nodes(:, 2));
ymax = max(mesh.Nodes(:, 2));
corners = [xmin ymin; xmax ymin; xmin ymax; xmax ymax];
labels = ["lower-left", "lower-right", "upper-left", "upper-right"];
[~, idx] = min(sum((corners - [x y]).^2, 2));
label = labels(idx);
end

function label = portable_source_label(studySet, caseName)
caseName = replace(string(caseName), "\\", "/");
if studySet == "MainCampaign"
    root = "Results_Parametric_ScientificallyCorrected_20260823";
elseif studySet == "PaperChecks"
    root = "Results_Paper_ScientificallyCorrected_20260823";
else
    root = string(studySet);
end
label = root + "/" + caseName + "/end.mat";
end

function write_readme(outputDir, summary, meshSummary)
reference = summary(summary.StudySet == "MainCampaign" & ...
    summary.Case == "Grid_G28x40", :);
fid = fopen(fullfile(outputDir, 'README_M7_STRAIN.md'), 'w');
if fid < 0
    error('M7Strain:ReadmeOpen', 'Could not create the M7 README.');
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# M7 endpoint-strain quantification\n\n');
fprintf(fid, 'Generated from saved `end.mat` states; no model case was rerun.\n\n');
fprintf(fid, '## Metric\n\n');
fprintf(fid, ['Infinitesimal von Mises-equivalent strain at all interior ' ...
    'Gauss points: `sqrt((2/3) dev(eps):dev(eps))`, with ' ...
    '`eps = sym(grad(u))` and plane-strain `eps_zz = 0`.\n\n']);
fprintf(fid, ['P99 is the whole-domain percentile. The raw maximum and its ' ...
    'coordinates are retained to expose constrained-corner localization. ' ...
    'The corner window is `x <= 0.20 mm, y <= 0.20 mm`.\n\n']);
fprintf(fid, ['Green-Lagrange comparisons are kinematic diagnostics only; ' ...
    'they are not a finite-strain validation of the coupled solution.\n\n']);
fprintf(fid, ['`det(F) <= 0` marks a kinematic inversion at an integration ' ...
    'point. Such raw maxima are reported for auditability but must not be ' ...
    'interpreted as physical material strains.\n\n']);
if ~isempty(reference)
    fprintf(fid, '## Reference result\n\n');
    fprintf(fid, ['For `Grid_G28x40`, whole-domain P99 = %.6g%% and the raw ' ...
        'maximum = %.6g%% at (%.6g, %.6g) mm.\n\n'], ...
        reference.EqvP99_pct, reference.EqvMax_pct, ...
        reference.EqvMaxX_mm, reference.EqvMaxY_mm);
end
fprintf(fid, '## Files\n\n');
fprintf(fid, '- `m7_strain_summary.csv`: all endpoint cases.\n');
fprintf(fid, '- `m7_mesh_strain_summary.csv`: three main-study meshes.\n');
fprintf(fid, ['- `m7_reference_Grid_G28x40_ip_strain.csv`: auditable ' ...
    'integration-point values for the reference grid.\n']);
fprintf(fid, '- `m7_strain_quantification.mat`: MATLAB tables.\n\n');
fprintf(fid, 'Three mesh rows written: %d.\n', height(meshSummary));
end
