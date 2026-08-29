classdef AgingDegradation < BaseModel
    %AGINGDEGRADATION Effective aging-state model for elastomer seals.
    %
    % Models the aging process of rubber/elastomer seals through:
    %   - Arrhenius-type thermal aging kinetics
    %   - Crosslink density evolution (hardening/embrittlement)
    %   - Transport of the phenomenological aging-extent field through the bulk
    %   - Spatial aging-extent propagation governed by diffusion and reaction
    %
    % The aging variable alpha_a ∈ [0,1] represents the degree of aging:
    %   0 = pristine material
    %   1 = fully aged/degraded material
    %
    % Evolution equation:
    %   d(alpha_a)/dt = k_age * (1 - alpha_a)^n_age * exp(-Ea/(R*T))
    %
    % where k_age is the aging rate constant, n_age is the reaction order,
    % Ea is the activation energy, and T is the temperature.
    %
    % Effective aging-extent transport and reaction:
    %   d(alpha_a)/dt = D_age * Laplacian(alpha_a) + reaction_term
    %
    % Required inputs:
    %   physics_in{k}.type      = "AgingDegradation";
    %   physics_in{k}.Egroup    = "Internal";
    %   physics_in{k}.k_age     = 1e-8;    % Aging rate constant [1/s]
    %   physics_in{k}.n_age     = 1.0;     % Reaction order [-]
    %   physics_in{k}.Ea_age    = 80e3;    % Activation energy [J/mol]
    %   physics_in{k}.D_age     = 1e-11;   % Effective aging-extent transport coefficient [m^2/s]
    %   physics_in{k}.T_service = 353.15;  % Service temperature [K] (80°C)
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        k_age       % Aging rate constant [1/s]
        n_age       % Reaction order [-]
        Ea_age      % Activation energy [J/mol]
        D_age       % Effective aging-extent transport coefficient [m^2/s]
        phiDiffusionFactor % Optional phase-field multiplier for aging transport [-]
        T_service   % Service temperature [K]
        
        R_const = 8.31446261815324;  % Gas constant
        T_ref = 293.15;              % Reference temperature
        
        aging_Step  % Solution step for aging variable
        phi_Step    % Phase-field step (for damage coupling)
        
        alpha_max   % Maximum aging level tracked
        alpha_avg   % Average aging level tracked
    end
    
    methods
        function obj = AgingDegradation(mesh, physics, inputs)
            obj.myName = "AgingDegradation";
            disp("Initializing " + obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"alpha_a", "phi"});
            obj.aging_Step = stp(1);
            obj.phi_Step = stp(2);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            obj.k_age = inputs.k_age;
            obj.n_age = inputs.n_age;
            obj.Ea_age = inputs.Ea_age;
            obj.D_age = inputs.D_age;
            if isfield(inputs, 'phiDiffusionFactor')
                obj.phiDiffusionFactor = inputs.phiDiffusionFactor;
            else
                obj.phiDiffusionFactor = 0.0;
            end
            obj.T_service = inputs.T_service;
            
            obj.alpha_max = 0;
            obj.alpha_avg = 0;
        end
        
        function getKf(obj, physics, stp)
            if stp ~= obj.aging_Step
                return;
            end
            
            fprintf("        AgingDegradation get Matrix:")
            t = tic;
            
            dt = physics.dt;
            SVec = physics.StateVec;
            SVecOld = physics.StateVec_Old;
            
            dofmatX = []; dofmatY = []; kmat = [];
            fvec = []; dofvec = [];
            
            % Temperature-dependent aging rate
            k_eff = obj.k_age * exp(-obj.Ea_age/obj.R_const * ...
                (1/obj.T_service - 1/obj.T_ref));
            
            alpha_sum = 0;
            alpha_max_local = 0;
            area_total = 0;
            
            for n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
                
                dofsAging = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                
                ALPHA = SVec{obj.aging_Step}(dofsAging);
                ALPHAOld = SVecOld{obj.aging_Step}(dofsAging);
                PHI = SVec{obj.phi_Step}(dofsPhi);
                
                f_el = zeros(length(dofsAging), 1);
                K_el = zeros(length(dofsAging));
                
                for ip=1:length(w)
                    alpha_ip = max(0, min(1, N(ip,:)*ALPHA));
                    alpha_old_ip = max(0, N(ip,:)*ALPHAOld);
                    phi_ip = min(max(N(ip,:)*PHI, 0), 1);
                    
                    %% Capacity: d(alpha)/dt
                    f_el = f_el + w(ip) * N(ip,:)'*N(ip,:)*(ALPHA - ALPHAOld)/dt;
                    K_el = K_el + w(ip) * N(ip,:)'*N(ip,:)/dt;
                    
                    %% Effective aging-extent transport
                    % The submission model uses a constant aging-extent
                    % diffusivity (phiDiffusionFactor = 0).  A nonzero multiplier
                    % is retained only as an explicit sensitivity option; the old
                    % hard-coded value of 1000 caused high-phi cases to be driven
                    % artificially toward the imposed alpha boundary value.
                    D_eff = obj.D_age * (1 + obj.phiDiffusionFactor*phi_ip);
                    f_el = f_el + w(ip) * D_eff * squeeze(G(ip,:,:))*squeeze(G(ip,:,:))'*ALPHA;
                    K_el = K_el + w(ip) * D_eff * squeeze(G(ip,:,:))*squeeze(G(ip,:,:))';
                    
                    %% Reaction term: k_eff * (1-alpha)^n
                    % Aging reaction (source term)
                    reaction = k_eff * max(0, (1 - alpha_ip))^obj.n_age;
                    
                    % Linearized: d(reaction)/d(alpha) for tangent
                    if alpha_ip < 1
                        d_reaction = -k_eff * obj.n_age * max(0, (1-alpha_ip))^(obj.n_age-1);
                    else
                        d_reaction = 0;
                    end
                    
                    f_el = f_el - w(ip) * N(ip,:)' * reaction;
                    K_el = K_el - w(ip) * N(ip,:)' * d_reaction * N(ip,:);
                    
                    % Track statistics
                    alpha_sum = alpha_sum + w(ip) * alpha_ip;
                    alpha_max_local = max(alpha_max_local, alpha_ip);
                    area_total = area_total + w(ip);
                end
                
                [dofmatxloc,dofmatyloc] = ndgrid(dofsAging, dofsAging);
                dofmatX = [dofmatX; dofmatxloc(:)];
                dofmatY = [dofmatY; dofmatyloc(:)];
                kmat = [kmat; K_el(:)];
                fvec = [fvec; f_el];
                dofvec = [dofvec; dofsAging];
            end
            
            physics.fint{obj.aging_Step} = physics.fint{obj.aging_Step} + ...
                sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{obj.aging_Step}), 1);
            physics.K{obj.aging_Step} = physics.K{obj.aging_Step} + ...
                sparse(dofmatX, dofmatY, kmat, length(physics.fint{obj.aging_Step}), length(physics.fint{obj.aging_Step}));
            
            if area_total > 0
                obj.alpha_avg = alpha_sum / area_total;
            end
            obj.alpha_max = alpha_max_local;
            
            tElapsed = toc(t);
            fprintf("            (Assemble time:" + string(tElapsed) + ")\n");
            fprintf("Aging: avg=%.4f, max=%.4f\n", obj.alpha_avg, obj.alpha_max);
        end
        
        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
            hasInfo = false;
            provided = [];
        end
    end
end
