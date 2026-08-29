classdef FluidDiffusion < BaseModel
    %FLUIDDIFFUSION Stable Fickian diffusion of hydraulic fluid into seal.
    %
    % Pure Fick's law: dC/dt = div(D_eff * grad(C))
    % Phase-field enhanced diffusivity in damaged regions.

    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices

        DL          % Bulk diffusivity [m^2/s]
        kmin        % Residual stiffness factor [-]
        DL_crack    % Enhanced diffusivity in damaged regions [m^2/s]

        CL_int      % Integrated fluid content
        CL_max      % Maximum fluid concentration

        dx_Step
        C_Step
        phi_Step
    end

    methods
        function obj = FluidDiffusion(mesh, physics, inputs)
            obj.myName = "FluidDiffusion";
            disp("Initializing " + obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;

            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"dx","dy","CL","phi"});
            obj.dx_Step = stp(1);
            obj.phi_Step = stp(4);
            obj.C_Step = stp(3);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));

            obj.DL = inputs.DL;
            obj.kmin = inputs.kmin;
            obj.DL_crack = inputs.DL_crack;
            % Note: C_sat, VF, Ea_diff, C10 are no longer used.
            % If passed in inputs struct, they are simply ignored.

            obj.CL_int = 0;
            obj.CL_max = 0;
        end

        function getKf(obj, physics, stp)
            if stp ~= obj.C_Step
                return;
            end

            fprintf("        FluidDiffusion get Matrix:")
            t = tic;

            dt = physics.dt;
            SVec = physics.StateVec;
            SVecOld = physics.StateVec_Old;

            dofmatX = []; dofmatY = []; kmat = [];
            fvec = []; dofvec = [];

            CL_sum = 0;
            CL_max2 = 0;

            for n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);

                dofsPHI = obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), Elem_Nodes);
                dofsCL = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);

                PHI = SVec{obj.phi_Step}(dofsPHI);
                CL = SVec{obj.C_Step}(dofsCL);
                CLOld = SVecOld{obj.C_Step}(dofsCL);

                q_el = zeros(length(dofsCL), 1);
                K_cc = zeros(length(dofsCL));

                for ip=1:length(w)
                    %% Capacity: dC/dt
                    q_el = q_el + w(ip) * N(ip,:)'*N(ip,:)*(CL-CLOld)/dt;
                    K_cc = K_cc + w(ip) * N(ip,:)'*N(ip,:)/dt;

                    %% Phase-field enhanced diffusivity
                    ff = min(max(N(ip,:)*PHI, 0), 1);
                    D_local = obj.DL*(1-ff) + obj.DL_crack*ff;

                    %% Pure Fickian diffusion: D * grad(C) . grad(N)
                    GG = squeeze(G(ip,:,:));
                    q_el = q_el + w(ip)*D_local * GG*GG'*CL;
                    K_cc = K_cc + w(ip)*D_local * GG*GG';

                    CL_sum = CL_sum + w(ip)*max(0, N(ip,:)*CL);
                    CL_max2 = max(CL_max2, N(ip,:)*CL);
                end

                [dofmatxloc,dofmatyloc] = ndgrid(dofsCL,dofsCL);
                dofmatX = [dofmatX; dofmatxloc(:)];
                dofmatY = [dofmatY; dofmatyloc(:)];
                kmat = [kmat; K_cc(:)];

                fvec = [fvec; q_el];
                dofvec = [dofvec; dofsCL];
            end

            physics.fint{obj.C_Step} = physics.fint{obj.C_Step} + ...
                sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{obj.C_Step}), 1);
            physics.K{obj.C_Step} = physics.K{obj.C_Step} + ...
                sparse(dofmatX, dofmatY, kmat, length(physics.fint{obj.C_Step}), length(physics.fint{obj.C_Step}));

            obj.CL_int = CL_sum;
            obj.CL_max = CL_max2;

            tElapsed = toc(t);
            fprintf("            (Assemble time:" + string(tElapsed) + ")\n");
        end

        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
            hasInfo = false;
            provided = [];
        end
    end
end
