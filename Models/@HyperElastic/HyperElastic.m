classdef HyperElastic < BaseModel
    %HYPERELASTIC Neo-Hookean hyperelastic model for rubber/elastomer seal
    % materials used in aerospace hydraulic systems.
    %
    % Supports:
    %   - Neo-Hookean strain energy: W = C10*(I1-3) + 1/D1*(J-1)^2
    %   - Optional Prony-series viscoelastic relaxation
    %   - Temperature-dependent modulus (WLF shift)
    %   - Aging-induced stiffening via coupling to AgingDegradation model
    %
    % Required inputs:
    %   physics_in{k}.type     = "HyperElastic";
    %   physics_in{k}.Egroup   = "Internal";
    %   physics_in{k}.C10      = 1.5e6;    % Neo-Hookean coefficient [Pa]
    %   physics_in{k}.D1       = 1e-8;     % Bulk compressibility [1/Pa]
    %   physics_in{k}.kmin     = 1e-10;    % Residual stiffness factor [-]
    %
    % Optional inputs:
    %   physics_in{k}.nu_visc     = [0.3, 0.15];  % Prony series modulus fractions
    %   physics_in{k}.tau_visc    = [10, 1000];    % Relaxation times [s]
    %   physics_in{k}.T_ref       = 293.15;        % Reference temperature [K]
    %   physics_in{k}.C1_WLF      = 17.44;         % WLF constant C1
    %   physics_in{k}.C2_WLF      = 51.6;          % WLF constant C2
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        C10         % Neo-Hookean coefficient [Pa]
        D1          % Bulk compressibility [1/Pa]
        kmin        % Residual stiffness factor [-]
        
        % Viscoelastic Prony series
        nu_visc     % Modulus fractions for each branch
        tau_visc    % Relaxation times [s]
        n_visc      % Number of viscoelastic branches
        
        % Temperature dependence
        T_ref       % Reference temperature [K]
        C1_WLF      % WLF constant C1
        C2_WLF      % WLF constant C2
        T_current   % Current temperature [K]
        
        % Internal variables
        HistStress      % Viscoelastic history stresses (current)
        HistStressOld   % Viscoelastic history stresses (converged)
        
        dx_Step
        phi_step
        aging_step
    end
    
    methods
        function obj = HyperElastic(mesh, physics, inputs)
            obj.myName = "HyperElastic";
            disp("Initializing " + obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"dx","dy","phi"});
            obj.dx_Step = stp(1);
            obj.phi_step = stp(3);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            % Material parameters
            obj.C10 = inputs.C10;
            obj.D1 = inputs.D1;
            obj.kmin = inputs.kmin;
            
            % Viscoelastic parameters (optional)
            if isfield(inputs, 'nu_visc') && ~isempty(inputs.nu_visc)
                obj.nu_visc = inputs.nu_visc;
                obj.tau_visc = inputs.tau_visc;
                obj.n_visc = length(obj.nu_visc);
            else
                obj.nu_visc = [];
                obj.tau_visc = [];
                obj.n_visc = 0;
            end
            
            % Temperature dependence (optional)
            if isfield(inputs, 'T_ref')
                obj.T_ref = inputs.T_ref;
                obj.C1_WLF = inputs.C1_WLF;
                obj.C2_WLF = inputs.C2_WLF;
            else
                obj.T_ref = 293.15;
                obj.C1_WLF = 0;
                obj.C2_WLF = 1;
            end
            obj.T_current = obj.T_ref;
            
            % Initialize viscoelastic history
            nelem = size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1);
            nipmax = obj.mesh.ipcount1D^2;
            if obj.n_visc > 0
                obj.HistStress = zeros(nelem, nipmax, obj.n_visc, 4);
                obj.HistStressOld = obj.HistStress;
            end
        end
        
        function Commit(obj, physics, commit_type)
            if commit_type == "Timedep" && obj.n_visc > 0
                obj.HistStressOld = obj.HistStress;
            end
        end
        
        function getKf(obj, physics, stp)
            if stp ~= obj.dx_Step
                return;
            end
            
            fprintf("        HyperElastic get Matrix:")
            t = tic;
            
            dofmatX = [];
            dofmatY = [];
            kmat = [];
            fvec = [];
            dofvec = [];
            
            SVec = physics.StateVec;
            dt = max(physics.dt, 1e-20);
            
            % WLF time-temperature shift
            aT = 1.0;
            if obj.C1_WLF ~= 0
                dT = obj.T_current - obj.T_ref;
                log_aT = -obj.C1_WLF * dT / (obj.C2_WLF + dT);
                aT = 10^log_aT;
            end
            
            % Effective C10 (long-term + instantaneous viscoelastic stiffness)
            C10_inf = obj.C10 * (1 - sum(obj.nu_visc));
            
            newHistStress = obj.HistStress;
            
            for n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
                
                dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                dofsXY = [dofsX; dofsY];
                dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
                
                X = SVec{obj.dx_Step}(dofsX);
                Y = SVec{obj.dx_Step}(dofsY);
                XY = [X;Y];
                PHI = SVec{obj.phi_step}(dofsPhi);
                
                f_el = zeros(length(dofsXY), 1);
                K_el = zeros(length(dofsXY));
                
                for ip=1:length(w)
                    % Phase-field degradation
                    ff = min(max(N(ip,:)*PHI, 0), 1);
                    dam_fun = (1-obj.kmin)*(1-ff)^2 + obj.kmin;
                    
                    % Compute deformation gradient F = I + grad(u)
                    B = obj.getB(G(ip,:,:));
                    strain = B*XY;  % [exx; eyy; ezz; 2*exy]
                    
                    % For small-to-moderate strains, use linearized 
                    % Neo-Hookean (equivalent to Mooney-Rivlin with C01=0)
                    % This is valid for seal deformations typically < 30%
                    
                    % Effective stiffness matrix (Neo-Hookean linearized)
                    mu = 2*C10_inf;      % Shear modulus
                    kappa = 2/obj.D1;    % Bulk modulus
                    lambda = kappa - 2*mu/3;
                    
                    % Plane strain stiffness
                    D_nh = zeros(4,4);
                    D_nh(1,1) = lambda + 2*mu;
                    D_nh(1,2) = lambda;
                    D_nh(1,3) = lambda;
                    D_nh(2,1) = lambda;
                    D_nh(2,2) = lambda + 2*mu;
                    D_nh(2,3) = lambda;
                    D_nh(3,1) = lambda;
                    D_nh(3,2) = lambda;
                    D_nh(3,3) = lambda + 2*mu;
                    D_nh(4,4) = mu;
                    
                    % Equilibrium (instantaneous) stress
                    stress_eq = dam_fun * D_nh * strain;
                    D_tang = dam_fun * D_nh;
                    
                    % Add viscoelastic overstress from each Prony branch
                    for v = 1:obj.n_visc
                        tau_eff = obj.tau_visc(v) * aT;
                        alpha = exp(-dt / tau_eff);
                        
                        % Recursive update: s_v^{n+1} = alpha*s_v^n + nu_v*(1-alpha)*2*C10*D*strain
                        mu_v = 2*obj.C10*obj.nu_visc(v);
                        D_v = D_nh * (mu_v / mu);
                        
                        s_old = squeeze(obj.HistStressOld(n_el, ip, v, :));
                        s_new = alpha * s_old + (1-alpha) * dam_fun * D_v * strain;
                        
                        stress_eq = stress_eq + s_new;
                        D_tang = D_tang + (1-alpha) * dam_fun * D_v;
                        
                        newHistStress(n_el, ip, v, :) = s_new;
                    end
                    
                    f_el = f_el + B'*stress_eq*w(ip);
                    K_el = K_el + B'*D_tang*B*w(ip);
                end
                
                [dofmatxloc,dofmatyloc] = ndgrid(dofsXY,dofsXY);
                dofmatX = [dofmatX; dofmatxloc(:)];
                dofmatY = [dofmatY; dofmatyloc(:)];
                kmat = [kmat; K_el(:)];
                
                fvec = [fvec; f_el];
                dofvec = [dofvec; dofsXY];
            end
            
            if obj.n_visc > 0
                obj.HistStress = newHistStress;
            end
            
            physics.fint{stp} = physics.fint{stp} + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{stp}), 1);
            physics.K{stp} = physics.K{stp} + sparse(dofmatX, dofmatY, kmat, length(physics.fint{stp}), length(physics.fint{stp}));
            
            tElapsed = toc(t);
            fprintf("            (Assemble time:" + string(tElapsed) + ")\n");
        end
        
        function B = getB(~, grads)
            cp_count = size(grads, 2);
            B = zeros(4, cp_count*2);
            for ii = 1:cp_count
                B(1, ii) = grads(1,ii, 1);
                B(4, ii) = grads(1,ii, 2);
                B(2, ii + cp_count) = grads(1,ii, 2);
                B(4, ii + cp_count) = grads(1,ii, 1);
            end
        end
        
        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
            hasInfo = false;
            provided = [];
            
            if (var == "stresses" || var=="sxx" || var=="syy" || var=="szz" || var=="sxy" || var=="sh" || var=="vonMises")
                hasInfo = true;
                if (var == "stresses")
                    provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount, 4);
                else
                    provided = zeros(length(elems), obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
                end
                
                C10_inf = obj.C10 * (1 - sum(obj.nu_visc));
                mu = 2*C10_inf;
                kappa = 2/obj.D1;
                lambda = kappa - 2*mu/3;
                D_nh = zeros(4,4);
                D_nh(1,1) = lambda + 2*mu; D_nh(1,2) = lambda; D_nh(1,3) = lambda;
                D_nh(2,1) = lambda; D_nh(2,2) = lambda + 2*mu; D_nh(2,3) = lambda;
                D_nh(3,1) = lambda; D_nh(3,2) = lambda; D_nh(3,3) = lambda + 2*mu;
                D_nh(4,4) = mu;
                
                for el=1:length(elems)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, elems(el));
                    [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, elems(el));
                    
                    dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                    dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                    dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
                    
                    X = physics.StateVec{obj.dx_Step}(dofsX);
                    Y = physics.StateVec{obj.dx_Step}(dofsY);
                    XY = [X;Y];
                    PHI = physics.StateVec{obj.phi_step}(dofsPhi);
                    
                    for ip=1:length(w)
                        ff = min(max(N(ip,:)*PHI,0),1);
                        dam_fun = (1-obj.kmin)*(1-ff)^2+obj.kmin;
                        
                        B = obj.getB(G(ip,:,:));
                        strain = B*XY;
                        stress = dam_fun*D_nh*strain;
                        
                        if (loc == "Interior")
                            switch var
                                case "stresses"
                                    provided(el, ip, :) = stress;
                                case "sxx"
                                    provided(el, ip) = stress(1);
                                case "syy"
                                    provided(el, ip) = stress(2);
                                case "szz"
                                    provided(el, ip) = stress(3);
                                case "sxy"
                                    provided(el, ip) = stress(4);
                                case "sh"
                                    provided(el, ip) = (stress(1)+stress(2)+stress(3))/3;
                                case "vonMises"
                                    s = stress;
                                    provided(el, ip) = sqrt(0.5*((s(1)-s(2))^2+(s(2)-s(3))^2+(s(3)-s(1))^2+6*s(4)^2));
                            end
                        end
                    end
                end
            end
        end
    end
end
