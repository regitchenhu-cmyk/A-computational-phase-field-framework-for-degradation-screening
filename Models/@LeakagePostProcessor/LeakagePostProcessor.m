classdef LeakagePostProcessor < BaseModel
    %LEAKAGEPOSTPROCESSOR Darcy leakage estimator driven by phase-field damage.
    %
    % This model is intentionally post-processing only: it reads the current
    % phase-field and a contact-closure proxy, solves a 2D Darcy pressure
    % problem from the high-pressure side to the low-pressure side, and
    % reports the leakage rate Q.
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        phiType
        phiStep
        
        P_high
        P_low
        mu
        k0
        k_phi
        phiExponent
        aperture0
        aperturePhi
        apertureExponent
        contactChi
        contactDecayLength
        contactPressure0
        contactPressureRef
        outOfPlaneWidth
        connectedPhiThreshold
        connectedLeakinessThreshold
        
        Q
        Q_poiseuille
        maxPermeability
        meanPermeability
        maxAperture
        meanAperture
        pathCost
        pathMeanAperture
        pathMinAperture
        contactPressureMean
        contactClosureMean
        connectedPath
        phiField
        leakinessField
        pressureField
        permeabilityField
        apertureField
        pathMask
        contactPressureField
    end
    
    methods
        function obj = LeakagePostProcessor(mesh, physics, inputs)
            obj.myName = "LeakagePostProcessor";
            disp("Initializing " + obj.myName)
            obj.mesh = mesh;
            obj.myGroup = get_input(inputs, 'Egroup', "Internal");
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            [obj.phiType, stp] = obj.dofSpace.getDofType({"phi"});
            obj.phiStep = stp(1);
            
            obj.P_high = get_input(inputs, 'P_high', 21e6);
            obj.P_low = get_input(inputs, 'P_low', 0.101e6);
            obj.mu = get_input(inputs, 'mu', 0.046);                  % calibrated hydraulic-fluid viscosity [Pa*s]
            obj.k0 = get_input(inputs, 'k0', 1e-22);                  % intact seal permeability [m^2]
            obj.k_phi = get_input(inputs, 'k_phi', 1e-14);            % damaged-path permeability scale [m^2]
            obj.phiExponent = get_input(inputs, 'phiExponent', 4.0);
            obj.aperture0 = get_input(inputs, 'aperture0', 1e-9);      % residual closed-channel aperture [m]
            obj.aperturePhi = get_input(inputs, 'aperturePhi', 5e-6);  % fully damaged aperture scale [m]
            obj.apertureExponent = get_input(inputs, 'apertureExponent', 2.0);
            obj.contactChi = get_input(inputs, 'contactChi', 3.0);
            obj.contactDecayLength = get_input(inputs, 'contactDecayLength', 0.5e-3);
            obj.contactPressure0 = get_input(inputs, 'contactPressure0', 1.2*obj.P_high);
            obj.contactPressureRef = get_input(inputs, 'contactPressureRef', obj.P_high);
            obj.outOfPlaneWidth = get_input(inputs, 'outOfPlaneWidth', 1.0);
            obj.connectedPhiThreshold = get_input(inputs, 'connectedPhiThreshold', 0.35);
            obj.connectedLeakinessThreshold = get_input(inputs, 'connectedLeakinessThreshold', 1e-3);
            
            obj.Q = 0;
            obj.Q_poiseuille = 0;
            obj.maxPermeability = obj.k0;
            obj.meanPermeability = obj.k0;
            obj.maxAperture = obj.aperture0;
            obj.meanAperture = obj.aperture0;
            obj.pathCost = inf;
            obj.pathMeanAperture = obj.aperture0;
            obj.pathMinAperture = obj.aperture0;
            obj.contactPressureMean = obj.contactPressure0;
            obj.contactClosureMean = 1;
            obj.connectedPath = false;
            obj.phiField = [];
            obj.leakinessField = [];
            obj.pressureField = [];
            obj.permeabilityField = [];
            obj.apertureField = [];
            obj.pathMask = [];
            obj.contactPressureField = [];
        end
        
        function metrics = Evaluate(obj, physics)
            allNodes = obj.mesh.GetAllNodesForGroup(obj.myGroupIndex);
            phiDofs = obj.dofSpace.getDofIndices(obj.phiType, allNodes);
            phiNodes = min(max(physics.StateVec{obj.phiStep}(phiDofs), 0), 1);
            
            coords = obj.mesh.Nodes(allNodes, :);
            [xVals, ~, ix] = unique(coords(:,1));
            [yVals, ~, iy] = unique(coords(:,2));
            nx = numel(xVals);
            ny = numel(yVals);
            
            phiGrid = zeros(ny, nx);
            for n = 1:numel(allNodes)
                phiGrid(iy(n), ix(n)) = max(phiGrid(iy(n), ix(n)), phiNodes(n));
            end
            
            closure = max(0, (1 - phiGrid).^2);
            pcBase = obj.estimateContactPressureGrid(physics, xVals, yVals);
            pc = pcBase .* closure;
            leakiness = phiGrid.^obj.phiExponent .* exp(-obj.contactChi * pc / max(obj.contactPressureRef, eps));
            kEff = obj.k0 + obj.k_phi * leakiness;
            aperture = obj.aperture0 + obj.aperturePhi * ...
                phiGrid.^obj.apertureExponent .* exp(-obj.contactChi * pc / max(obj.contactPressureRef, eps));
            
            [Qval, pGrid] = obj.solveDarcy(xVals, yVals, kEff);
            Qp = obj.solvePoiseuille(xVals, yVals, aperture);
            path = obj.findMinimumLeakagePath(kEff, aperture);
            
            obj.Q = Qval;
            obj.Q_poiseuille = Qp;
            obj.maxPermeability = max(kEff(:));
            obj.meanPermeability = mean(kEff(:));
            obj.maxAperture = max(aperture(:));
            obj.meanAperture = mean(aperture(:));
            obj.pathCost = path.cost;
            obj.pathMeanAperture = path.meanAperture;
            obj.pathMinAperture = path.minAperture;
            obj.contactPressureMean = mean(pc(:));
            obj.contactClosureMean = mean(closure(:));
            obj.connectedPath = obj.hasConnectedPath(phiGrid, leakiness);
            obj.phiField = phiGrid;
            obj.leakinessField = leakiness;
            obj.pressureField = pGrid;
            obj.permeabilityField = kEff;
            obj.apertureField = aperture;
            obj.pathMask = path.mask;
            obj.contactPressureField = pc;
            
            metrics.Q = obj.Q;
            metrics.Q_poiseuille = obj.Q_poiseuille;
            metrics.maxPermeability = obj.maxPermeability;
            metrics.meanPermeability = obj.meanPermeability;
            metrics.maxAperture = obj.maxAperture;
            metrics.meanAperture = obj.meanAperture;
            metrics.pathCost = obj.pathCost;
            metrics.pathMeanAperture = obj.pathMeanAperture;
            metrics.pathMinAperture = obj.pathMinAperture;
            metrics.contactPressureMean = obj.contactPressureMean;
            metrics.contactClosureMean = obj.contactClosureMean;
            metrics.connectedPath = obj.connectedPath;
        end
        
        function [Qval, pGrid] = solveDarcy(obj, xVals, yVals, kEff)
            nx = numel(xVals);
            ny = numel(yVals);
            nUnknown = nx * ny;
            idx = @(j,i) (i-1)*ny + j;
            A = spalloc(nUnknown, nUnknown, 5*nUnknown);
            b = zeros(nUnknown, 1);
            kSolve = kEff / max(max(kEff(:)), eps);
            % Omitting neighbors outside j=1 and j=ny implements zero
            % normal Darcy flux on the lower and upper edges.
            
            for i = 1:nx
                for j = 1:ny
                    row = idx(j,i);
                    if i == 1
                        A(row,row) = 1;
                        b(row) = obj.P_high;
                        continue;
                    elseif i == nx
                        A(row,row) = 1;
                        b(row) = obj.P_low;
                        continue;
                    end
                    
                    diagVal = 0;
                    if i > 1
                        g = obj.edgeConductance(kSolve(j,i), kSolve(j,i-1), xVals(i)-xVals(i-1), obj.localSpacing(yVals,j));
                        A(row, idx(j,i-1)) = A(row, idx(j,i-1)) - g;
                        diagVal = diagVal + g;
                    end
                    if i < nx
                        g = obj.edgeConductance(kSolve(j,i), kSolve(j,i+1), xVals(i+1)-xVals(i), obj.localSpacing(yVals,j));
                        A(row, idx(j,i+1)) = A(row, idx(j,i+1)) - g;
                        diagVal = diagVal + g;
                    end
                    if j > 1
                        g = obj.edgeConductance(kSolve(j,i), kSolve(j-1,i), yVals(j)-yVals(j-1), obj.localSpacing(xVals,i));
                        A(row, idx(j-1,i)) = A(row, idx(j-1,i)) - g;
                        diagVal = diagVal + g;
                    end
                    if j < ny
                        g = obj.edgeConductance(kSolve(j,i), kSolve(j+1,i), yVals(j+1)-yVals(j), obj.localSpacing(xVals,i));
                        A(row, idx(j+1,i)) = A(row, idx(j+1,i)) - g;
                        diagVal = diagVal + g;
                    end
                    A(row,row) = diagVal;
                end
            end
            
            p = A \ b;
            pGrid = reshape(p, ny, nx);
            
            Qval = 0;
            dx = xVals(nx) - xVals(nx-1);
            for j = 1:ny
                dyw = obj.localSpacing(yVals, j);
                kFace = obj.harmonicMean(kEff(j,nx), kEff(j,nx-1));
                qj = kFace / obj.mu * (pGrid(j,nx-1) - pGrid(j,nx)) / dx * dyw * obj.outOfPlaneWidth;
                Qval = Qval + max(qj, 0);
            end
        end
        
        function Qval = solvePoiseuille(obj, xVals, yVals, aperture)
            L = max(xVals) - min(xVals);
            dp = max(obj.P_high - obj.P_low, 0);
            Qval = 0;
            if L <= 0 || dp <= 0
                return;
            end
            
            for j = 1:numel(yVals)
                dyw = obj.localSpacing(yVals, j);
                hRow = max(aperture(j, :), [], 2);
                % dyw already supplies the slit-width measure for the row
                % integral; multiplying by outOfPlaneWidth again would give
                % Q_poiseuille dimensions of m^4/s instead of m^3/s.
                Qval = Qval + hRow^3 / (12 * obj.mu) * dp / L * dyw;
            end
        end
        
        function g = edgeConductance(obj, k1, k2, distance, areaLength)
            kFace = obj.harmonicMean(k1, k2);
            g = kFace / obj.mu * areaLength * obj.outOfPlaneWidth / max(distance, eps);
        end
        
        function h = harmonicMean(~, a, b)
            h = 2*a*b / max(a+b, eps);
        end
        
        function ds = localSpacing(~, vals, i)
            n = numel(vals);
            if n == 1
                ds = 1;
            elseif i == 1
                ds = 0.5 * (vals(2) - vals(1));
            elseif i == n
                ds = 0.5 * (vals(n) - vals(n-1));
            else
                ds = 0.5 * (vals(i+1) - vals(i-1));
            end
            ds = max(ds, eps);
        end
        
        function isConnected = hasConnectedPath(obj, phiGrid, leakiness)
            open = phiGrid >= obj.connectedPhiThreshold | leakiness >= obj.connectedLeakinessThreshold;
            [ny, nx] = size(open);
            visited = false(ny, nx);
            queue = zeros(numel(open), 2);
            head = 1;
            tail = 0;
            
            for j = 1:ny
                if open(j,1)
                    tail = tail + 1;
                    queue(tail,:) = [j, 1];
                    visited(j,1) = true;
                end
            end
            
            isConnected = false;
            while head <= tail
                j = queue(head,1);
                i = queue(head,2);
                head = head + 1;
                if i == nx
                    isConnected = true;
                    return;
                end
                neigh = [j-1 i; j+1 i; j i-1; j i+1];
                for n = 1:4
                    jj = neigh(n,1);
                    ii = neigh(n,2);
                    if jj >= 1 && jj <= ny && ii >= 1 && ii <= nx && open(jj,ii) && ~visited(jj,ii)
                        tail = tail + 1;
                        queue(tail,:) = [jj, ii];
                        visited(jj,ii) = true;
                    end
                end
            end
        end
        
        function path = findMinimumLeakagePath(obj, kEff, aperture)
            [ny, nx] = size(kEff);
            conductance = kEff / max(max(kEff(:)), eps);
            nodeCost = 1 ./ max(conductance, eps);
            dist = inf(ny, nx);
            prevJ = zeros(ny, nx);
            prevI = zeros(ny, nx);
            visited = false(ny, nx);
            
            dist(:,1) = nodeCost(:,1);
            for iter = 1:numel(kEff)
                candidate = dist;
                candidate(visited) = inf;
                [bestCost, linearIdx] = min(candidate(:));
                if ~isfinite(bestCost)
                    break;
                end
                [j, i] = ind2sub([ny, nx], linearIdx);
                visited(j,i) = true;
                if i == nx
                    break;
                end
                
                neigh = [j-1 i; j+1 i; j i-1; j i+1];
                for n = 1:4
                    jj = neigh(n,1);
                    ii = neigh(n,2);
                    if jj < 1 || jj > ny || ii < 1 || ii > nx || visited(jj,ii)
                        continue;
                    end
                    stepCost = 0.5 * (nodeCost(j,i) + nodeCost(jj,ii));
                    newCost = dist(j,i) + stepCost;
                    if newCost < dist(jj,ii)
                        dist(jj,ii) = newCost;
                        prevJ(jj,ii) = j;
                        prevI(jj,ii) = i;
                    end
                end
            end
            
            [pathCost, endJ] = min(dist(:,nx));
            endI = nx;
            mask = false(ny, nx);
            apertureVals = [];
            
            if isfinite(pathCost)
                j = endJ;
                i = endI;
                while i >= 1 && j >= 1
                    mask(j,i) = true;
                    apertureVals(end+1) = aperture(j,i); %#ok<AGROW>
                    if i == 1
                        break;
                    end
                    pj = prevJ(j,i);
                    pi = prevI(j,i);
                    if pj == 0 || pi == 0
                        break;
                    end
                    j = pj;
                    i = pi;
                end
            end
            
            path.cost = pathCost;
            path.mask = mask;
            if isempty(apertureVals)
                path.meanAperture = 0;
                path.minAperture = 0;
            else
                path.meanAperture = mean(apertureVals);
                path.minAperture = min(apertureVals);
            end
        end
        
        function pcGrid = estimateContactPressureGrid(obj, physics, xVals, yVals)
            nx = numel(xVals);
            ny = numel(yVals);
            [X, Y] = meshgrid(xVals, yVals);
            xmin = min(xVals);
            xmax = max(xVals);
            ymin = min(yVals);
            ymax = max(yVals);
            decay = max(obj.contactDecayLength, eps);
            
            pcGrid = zeros(ny, nx);
            foundContact = false;
            
            for m = 1:length(physics.models)
                mdl = physics.models{m};
                if ~(isprop(mdl, 'myName') && mdl.myName == "ContactPenalty" && ...
                        isprop(mdl, 'nodalPressure') && ~isempty(mdl.nodalPressure))
                    continue;
                end
                nodes = mdl.mesh.Nodegroups{mdl.myGroupIndex}.Nodes;
                coords = mdl.mesh.Nodes(nodes, :);
                pressure = mdl.nodalPressure(:);
                if isempty(pressure) || max(pressure) <= 0
                    continue;
                end
                foundContact = true;
                
                switch mdl.myGroup
                    case "Top"
                        pLine = obj.interpolateContactLine(coords(:,1), pressure, xVals);
                        pcGrid = max(pcGrid, repmat(pLine(:)', ny, 1) .* exp(-(ymax - Y) / decay));
                    case "Bottom"
                        pLine = obj.interpolateContactLine(coords(:,1), pressure, xVals);
                        pcGrid = max(pcGrid, repmat(pLine(:)', ny, 1) .* exp(-(Y - ymin) / decay));
                    case "Right"
                        pLine = obj.interpolateContactLine(coords(:,2), pressure, yVals);
                        pcGrid = max(pcGrid, repmat(pLine(:), 1, nx) .* exp(-(xmax - X) / decay));
                    case "Left"
                        pLine = obj.interpolateContactLine(coords(:,2), pressure, yVals);
                        pcGrid = max(pcGrid, repmat(pLine(:), 1, nx) .* exp(-(X - xmin) / decay));
                end
            end
            
            if ~foundContact
                pcGrid(:) = obj.contactPressure0;
            end
        end
        
        function pLine = interpolateContactLine(~, coord, pressure, queryCoord)
            [coordSorted, order] = sort(coord(:));
            pressureSorted = pressure(order);
            [coordUnique, ~, groupId] = unique(coordSorted);
            pressureUnique = accumarray(groupId, pressureSorted, [], @mean);
            
            if numel(coordUnique) == 1
                pLine = zeros(size(queryCoord)) + pressureUnique(1);
            else
                pLine = interp1(coordUnique, pressureUnique, queryCoord, 'linear', 'extrap');
            end
            pLine = max(pLine, 0);
        end
        
        function pcBase = estimateContactPressureBase(obj, physics)
            pcVals = [];
            for m = 1:length(physics.models)
                mdl = physics.models{m};
                if isprop(mdl, 'myName') && mdl.myName == "ContactPenalty" && ...
                        isprop(mdl, 'meanPressure') && mdl.meanPressure > 0
                    pcVals(end+1) = mdl.meanPressure; %#ok<AGROW>
                end
            end
            if isempty(pcVals)
                pcBase = obj.contactPressure0;
            else
                pcBase = max(pcVals);
            end
        end
    end
end

function value = get_input(s, name, defaultValue)
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end
