classdef PressureLoad < BaseModel
    %PRESSURELOAD Applies hydraulic pressure as a distributed load on
    % boundary elements. Supports time-varying (cyclic) pressure profiles.
    %
    % Required inputs:
    %   physics_in{k}.type       = "PressureLoad";
    %   physics_in{k}.Egroup     = "Left";       % Boundary element group
    %   physics_in{k}.pressure   = 21e6;         % Pressure magnitude [Pa]
    %   physics_in{k}.direction  = "normal";     % "normal" or [nx, ny]
    %
    % Optional inputs:
    %   physics_in{k}.P_profile  = "cyclic";     % "constant","cyclic","ramp"
    %   physics_in{k}.P_min      = 0;            % Min pressure for cyclic [Pa]
    %   physics_in{k}.f_cycle    = 0.1;          % Cycling frequency [Hz]
    %   physics_in{k}.t_ramp     = 10;           % Ramp time [s]
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        
        pressure        % Pressure magnitude [Pa]
        direction       % Load direction: "normal" or [nx, ny]
        P_profile       % Pressure profile type
        P_min           % Minimum pressure (cyclic) [Pa]
        f_cycle         % Cycling frequency [Hz]
        t_ramp          % Ramp time [s]
        
        dx_Step
    end
    
    methods
        function obj = PressureLoad(mesh, physics, inputs)
            obj.myName = "PressureLoad";
            disp("Initializing " + obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"dx","dy"});
            obj.dx_Step = stp(1);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            obj.pressure = inputs.pressure;
            obj.direction = inputs.direction;
            
            if isfield(inputs, 'P_profile')
                obj.P_profile = inputs.P_profile;
            else
                obj.P_profile = "constant";
            end
            
            if isfield(inputs, 'P_min')
                obj.P_min = inputs.P_min;
            else
                obj.P_min = 0;
            end
            
            if isfield(inputs, 'f_cycle')
                obj.f_cycle = inputs.f_cycle;
            else
                obj.f_cycle = 0.1;
            end
            
            if isfield(inputs, 't_ramp')
                obj.t_ramp = inputs.t_ramp;
            else
                obj.t_ramp = 10;
            end
        end
        
        function getKf(obj, physics, stp)
            if stp ~= obj.dx_Step
                return;
            end
            
            fprintf("        PressureLoad get force vector:")
            t_tic = tic;
            
            % Determine current pressure based on profile
            t_current = physics.time;
            switch obj.P_profile
                case "constant"
                    P_now = obj.pressure;
                case "ramp"
                    P_now = obj.pressure * min(1, t_current / obj.t_ramp);
                case "cyclic"
                    P_amp = (obj.pressure - obj.P_min) / 2;
                    P_mean = (obj.pressure + obj.P_min) / 2;
                    P_now = P_mean + P_amp * sin(2*pi*obj.f_cycle*t_current);
                otherwise
                    P_now = obj.pressure;
            end
            
            fvec = [];
            dofvec = [];
            
            for n_el = 1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
                [normals, ~] = obj.mesh.getNormals(obj.myGroupIndex, n_el);
                
                dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                
                fx_el = zeros(length(dofsX), 1);
                fy_el = zeros(length(dofsY), 1);
                
                for ip = 1:length(w)
                    if obj.direction == "normal"
                        % Outward normal (points out of the body)
                        nx = normals(ip, 1);
                        ny = normals(ip, 2);
                    else
                        nx = obj.direction(1);
                        ny = obj.direction(2);
                    end
                    
                    % Pressure traction on the boundary (Cauchy):
                    %     t = -p * n_outward    (compressive on the body)
                    % External nodal force:  f_ext = ∫ N · t dΓ = -p ∫ N n dΓ
                    % Solver residual is  R = f_int - f_ext = 0,  so we
                    % accumulate  -f_ext  into f_int:
                    %     f_int += -f_ext = +p ∫ N n dΓ
                    fx_el = fx_el + P_now * nx * N(ip,:)' * w(ip);
                    fy_el = fy_el + P_now * ny * N(ip,:)' * w(ip);
                end
                
                fvec = [fvec; fx_el; fy_el];
                dofvec = [dofvec; dofsX; dofsY];
            end
            
            % Pressure is an external load, so it enters as negative fint
            % (since solver solves: fint = 0, and fext enters as -fext in fint)
            physics.fint{stp} = physics.fint{stp} + ...
                sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{stp}), 1);
            
            tElapsed = toc(t_tic);
            fprintf("            (Assemble time:" + string(tElapsed) + ")\n");
            fprintf("Current pressure: %.2f MPa\n", P_now/1e6);
        end
        
        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
            hasInfo = false;
            provided = [];
        end
    end
end
