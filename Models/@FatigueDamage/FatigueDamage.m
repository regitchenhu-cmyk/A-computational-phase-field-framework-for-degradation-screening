classdef FatigueDamage < BaseModel
    %FATIGUEDAMAGE Capped cycle-based fatigue state for numerical screening.
    %
    % The local energetic branch is activated by the tensile-energy threshold,
    % combined with the pressure/contact-scaled base branch, and capped by maxRate.
    % With the reported parameters, every activated integration point uses the capped rate.
    % The committed integration-point state is projected to nodes by L2 projection.
    %
    % NOTE: this version DOES NOT solve a PDE for D_f. The damage
    % accumulates POINTWISE based on the local strain energy and is
    % then projected (least squares) onto the nodal dofs via a mass
    % matrix. This avoids the spurious diffusion-style smoothing that
    % was previously contaminating the solution and caused D_f > 1.
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        W0              % Fatigue threshold strain energy [J/m^3]
        beta            % Paris exponent
        f_cycle         % Cycling frequency [Hz]
        P_max
        P_min
        P_ref
        pressureExponent
        energyRateScale % Energetic fatigue-rate scale r0 [cycle^-1]
        contactFatigueCoeff
        contactFatigueExponent
        contactPressure0
        contactPressureRef
        contactDecayLength
        baseRate
        maxRate
        kmin
        
        Df_ip           % committed integration-point damage history [nelem x nip]
        Df_ip_trial     % current-step trial history; committed once after convergence
        Df_max          % maximum integration-point fatigue damage
        totalCycles
        
        dx_Step
        Df_Step
        phi_Step
    end
    
    methods
        function obj = FatigueDamage(mesh, physics, inputs)
            obj.myName = "FatigueDamage";
            disp("Initializing " + obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"dx","dy","D_f","phi"});
            obj.dx_Step = stp(1);
            obj.Df_Step = stp(3);
            obj.phi_Step = stp(4);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            obj.W0      = inputs.W0;
            obj.beta    = inputs.beta;
            obj.f_cycle = inputs.f_cycle;
            obj.P_max   = inputs.P_max;
            obj.P_min   = inputs.P_min;
            obj.P_ref   = get_input(inputs, 'P_ref', 21e6);
            obj.pressureExponent = get_input(inputs, 'pressureExponent', 0.0);
            obj.energyRateScale = get_input(inputs, 'energyRateScale', 1.0);
            obj.contactFatigueCoeff = get_input(inputs, 'contactFatigueCoeff', 0.0);
            obj.contactFatigueExponent = get_input(inputs, 'contactFatigueExponent', 1.0);
            obj.contactPressure0 = get_input(inputs, 'contactPressure0', obj.P_max);
            obj.contactPressureRef = get_input(inputs, 'contactPressureRef', obj.P_ref);
            obj.contactDecayLength = get_input(inputs, 'contactDecayLength', 0.5e-3);
            obj.baseRate = get_input(inputs, 'baseRate', 0.0);
            obj.maxRate = get_input(inputs, 'maxRate', 1e-5);
            obj.kmin    = inputs.kmin;
            
            nelem = size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1);
            nipmax = obj.mesh.ipcount1D^2;
            obj.Df_ip = zeros(nelem, nipmax);
            obj.Df_ip_trial = obj.Df_ip;
            obj.Df_max = 0;
            obj.totalCycles = 0;
        end
        
        function Commit(obj, physics, commit_type)
            if commit_type == "Timedep"
                % Assemble can be called repeatedly by Newton and staggered
                % iterations.  Only the final trial state is committed once
                % per converged physical time increment.
                if isempty(obj.Df_ip_trial) || ~isequal(size(obj.Df_ip_trial), size(obj.Df_ip))
                    obj.Df_ip_trial = obj.Df_ip;
                end
                obj.Df_ip = obj.Df_ip_trial;
                obj.Df_max = max(obj.Df_ip(:));
                dt = physics.dt;
                obj.totalCycles = obj.totalCycles + dt * obj.f_cycle;
            end
        end
        
        function getKf(obj, physics, stp)
            % Solve a simple L2 projection: (M) D_f_node = (M_lump) D_f_ip
            % This produces a smooth nodal field while keeping D_f bounded.
            if stp ~= obj.Df_Step
                return;
            end
            
            fprintf("        FatigueDamage:")
            t = tic;
            
            dt = physics.dt;
            SVec = physics.StateVec;
            dN_step = max(0, dt * obj.f_cycle);
            
            % Material moduli for tensile energy
            C10 = 1.5e6; D1 = 1e-9;
            for m = 1:length(physics.models)
                if isprop(physics.models{m}, 'C10') && ~isempty(physics.models{m}.C10)
                    C10 = physics.models{m}.C10;
                    D1  = physics.models{m}.D1;
                    break;
                end
            end
            mu     = 2*C10;
            kappa  = 2/D1;
            lambda = kappa - 2*mu/3;
            
            dofmatX = []; dofmatY = []; kmat = [];
            fvec = []; dofvec = [];
            
            nElem = size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1);
            newDf_ip = obj.Df_ip;
            contactPressureBase = obj.estimateContactPressureBase(physics);
            contactCache = obj.buildContactPressureCache(physics);
            
            for n_el=1:nElem
                Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
                
                dofsX  = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY  = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                dofsDf = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
                dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), Elem_Nodes);
                
                XY = [SVec{obj.dx_Step}(dofsX); SVec{obj.dx_Step}(dofsY)];
                PHI = SVec{obj.phi_Step}(dofsPhi);
                
                M_el  = zeros(length(dofsDf));
                rhs_el = zeros(length(dofsDf), 1);
                
                for ip=1:length(w)
                    B = getB_local(G(ip,:,:));
                    strain = B*XY;
                    xy_ip = N(ip,:) * obj.mesh.Nodes(Elem_Nodes, :);
                    
                    % --- Tensile strain energy (Miehe split) ---
                    exx=strain(1); eyy=strain(2); exy=strain(4)/2;
                    tr_e = exx+eyy;
                    e_avg=(exx+eyy)/2;
                    R = sqrt(((exx-eyy)/2)^2 + exy^2);
                    e1 = e_avg + R; e2 = e_avg - R;
                    Wplus = lambda/2*max(tr_e,0)^2 + mu*(max(e1,0)^2 + max(e2,0)^2);
                    
                    % --- Per-cycle damage rate (capped) ---
                    pressureAmp = max(0, obj.P_max - obj.P_min);
                    pressureFactor = max(0, pressureAmp / max(obj.P_ref, eps))^obj.pressureExponent;
                    baseTerm = obj.baseRate * max(0, pressureAmp / max(obj.P_ref, eps))^obj.beta;
                    if Wplus > 0.05*obj.W0
                        dDdN_ip = obj.energyRateScale * pressureFactor * ...
                            (Wplus / obj.W0)^obj.beta + baseTerm;
                    else
                        dDdN_ip = baseTerm;
                    end
                    phi_ip = min(max(N(ip,:) * PHI, 0), 1);
                    contactPressure = obj.estimateContactPressureAtPointCached(xy_ip, contactPressureBase, contactCache) * ...
                        max(0, (1 - phi_ip)^2);
                    contactFactor = 1 + obj.contactFatigueCoeff * ...
                        max(0, contactPressure / max(obj.contactPressureRef, eps))^obj.contactFatigueExponent;
                    dDdN_ip = dDdN_ip * contactFactor;
                    dDdN_ip = min(dDdN_ip, obj.maxRate); % per-cycle rate cap
                    
                    % --- Pointwise damage update ---
                    Df_old = obj.Df_ip(n_el, ip);
                    Df_new = min(Df_old + dDdN_ip * dN_step, 1.0);
                    newDf_ip(n_el, ip) = Df_new;
                    
                    % --- L2 projection: M * Df_node = N' * Df_ip ---
                    M_el  = M_el  + w(ip) * (N(ip,:)' * N(ip,:));
                    Df_new = double(Df_new(1));  % 强制标量
                    rhs_el = rhs_el + w(ip) * N(ip,:)' * Df_new;
                end
                
                % Equation we assemble: M * Df_node - rhs = 0
                f_el = M_el * SVec{obj.Df_Step}(dofsDf) - rhs_el;
                K_el = M_el;
                
                [dofmatxloc,dofmatyloc] = ndgrid(dofsDf,dofsDf);
                dofmatX = [dofmatX; dofmatxloc(:)];
                dofmatY = [dofmatY; dofmatyloc(:)];
                kmat    = [kmat; K_el(:)];
                fvec    = [fvec; f_el];
                dofvec  = [dofvec; dofsDf];
            end
            
            obj.Df_ip_trial = newDf_ip;
            Df_trial_max = max(newDf_ip(:));
            
            physics.fint{obj.Df_Step} = physics.fint{obj.Df_Step} + ...
                sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{obj.Df_Step}), 1);
            physics.K{obj.Df_Step} = physics.K{obj.Df_Step} + ...
                sparse(dofmatX, dofmatY, kmat, length(physics.fint{obj.Df_Step}), length(physics.fint{obj.Df_Step}));
            
            tElapsed = toc(t);
            fprintf("  (%.2fs)  Df_ip_max=%.4f  Cycles=%.0f\n", ...
                tElapsed, Df_trial_max, obj.totalCycles);
        end
        
        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
            hasInfo = false;
            provided = [];
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
        
        function pc = estimateContactPressureAtPoint(obj, physics, xy, fallbackPc)
            decay = max(obj.contactDecayLength, eps);
            pc = 0;
            foundContact = false;
            
            for m = 1:length(physics.models)
                mdl = physics.models{m};
                if ~(isprop(mdl, 'myName') && mdl.myName == "ContactPenalty" && ...
                        isprop(mdl, 'nodalPressure') && ~isempty(mdl.nodalPressure))
                    continue;
                end
                pressure = mdl.nodalPressure(:);
                if isempty(pressure) || max(pressure) <= 0
                    continue;
                end
                nodes = mdl.mesh.Nodegroups{mdl.myGroupIndex}.Nodes;
                coords = mdl.mesh.Nodes(nodes, :);
                foundContact = true;
                
                switch mdl.myGroup
                    case {"Top", "Bottom"}
                        pLine = interpolate_contact_line(coords(:,1), pressure, xy(1));
                        distance = abs(xy(2) - mean(coords(:,2)));
                    case {"Right", "Left"}
                        pLine = interpolate_contact_line(coords(:,2), pressure, xy(2));
                        distance = abs(xy(1) - mean(coords(:,1)));
                    otherwise
                        pLine = mean(pressure);
                        distance = 0;
                end
                pc = max(pc, max(pLine, 0) * exp(-distance / decay));
            end
            
            if ~foundContact
                pc = fallbackPc;
            end
        end
        
        function cache = buildContactPressureCache(obj, physics)
            cache = struct('group', {}, 'coord', {}, 'pressure', {}, 'boundaryCoord', {});
            for m = 1:length(physics.models)
                mdl = physics.models{m};
                if ~(isprop(mdl, 'myName') && mdl.myName == "ContactPenalty" && ...
                        isprop(mdl, 'nodalPressure') && ~isempty(mdl.nodalPressure))
                    continue;
                end
                pressure = mdl.nodalPressure(:);
                if isempty(pressure) || max(pressure) <= 0
                    continue;
                end
                nodes = mdl.mesh.Nodegroups{mdl.myGroupIndex}.Nodes;
                coords = mdl.mesh.Nodes(nodes, :);
                switch mdl.myGroup
                    case {"Top", "Bottom"}
                        lineCoord = coords(:,1);
                        boundaryCoord = mean(coords(:,2));
                    case {"Right", "Left"}
                        lineCoord = coords(:,2);
                        boundaryCoord = mean(coords(:,1));
                    otherwise
                        lineCoord = coords(:,1);
                        boundaryCoord = 0;
                end
                [coordUnique, pressureUnique] = prepare_contact_line(lineCoord, pressure);
                cache(end+1).group = mdl.myGroup; %#ok<AGROW>
                cache(end).coord = coordUnique;
                cache(end).pressure = pressureUnique;
                cache(end).boundaryCoord = boundaryCoord;
            end
        end
        
        function pc = estimateContactPressureAtPointCached(obj, xy, fallbackPc, cache)
            decay = max(obj.contactDecayLength, eps);
            pc = 0;
            for c = 1:numel(cache)
                switch cache(c).group
                    case {"Top", "Bottom"}
                        pLine = interp_contact_line(cache(c).coord, cache(c).pressure, xy(1));
                        distance = abs(xy(2) - cache(c).boundaryCoord);
                    case {"Right", "Left"}
                        pLine = interp_contact_line(cache(c).coord, cache(c).pressure, xy(2));
                        distance = abs(xy(1) - cache(c).boundaryCoord);
                    otherwise
                        pLine = mean(cache(c).pressure);
                        distance = 0;
                end
                pc = max(pc, max(pLine, 0) * exp(-distance / decay));
            end
            if isempty(cache)
                pc = fallbackPc;
            end
        end
    end
end

function value = get_input(inputs, name, defaultValue)
    if isfield(inputs, name) && ~isempty(inputs.(name))
        value = inputs.(name);
    else
        value = defaultValue;
    end
end

function B = getB_local(grads)
    cp_count = size(grads, 2);
    B = zeros(4, cp_count*2);
    for ii = 1:cp_count
        B(1, ii) = grads(1,ii, 1);
        B(4, ii) = grads(1,ii, 2);
        B(2, ii + cp_count) = grads(1,ii, 2);
        B(4, ii + cp_count) = grads(1,ii, 1);
    end
end

function p = interpolate_contact_line(coord, pressure, queryCoord)
    [coordUnique, pressureUnique] = prepare_contact_line(coord, pressure);
    p = interp_contact_line(coordUnique, pressureUnique, queryCoord);
end

function [coordUnique, pressureUnique] = prepare_contact_line(coord, pressure)
    [coordSorted, order] = sort(coord(:));
    pressureSorted = pressure(order);
    [coordUnique, ~, groupId] = unique(coordSorted);
    pressureUnique = accumarray(groupId, pressureSorted, [], @mean);
end

function p = interp_contact_line(coordUnique, pressureUnique, queryCoord)
    if numel(coordUnique) == 1
        p = pressureUnique(1);
    else
        p = interp1(coordUnique, pressureUnique, queryCoord, 'linear', 'extrap');
    end
    p = max(p, 0);
end
