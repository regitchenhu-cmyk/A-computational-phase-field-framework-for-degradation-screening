classdef ContactPenalty < BaseModel
    %CONTACTPENALTY Node-based rigid-wall contact for seal-groove studies.
    %
    % This is a lightweight contact-aware boundary model. It replaces hard
    % displacement constraints with a normal penalty reaction against a
    % rigid wall or platen and exposes mean/max contact pressure metrics.
    
    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices
        dx_Step
        
        axisName
        side
        wallPosition
        penalty
        pressurePenalty
        activeTolerance
        nodeArea
        
        meanPressure
        maxPressure
        activeFraction
        nodalPressure
        nodalPenetration
    end
    
    methods
        function obj = ContactPenalty(mesh, physics, inputs)
            obj.myName = "ContactPenalty";
            disp("Initializing " + obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Ngroup;
            obj.myGroupIndex = obj.mesh.getNodeGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"dx","dy"});
            obj.dx_Step = stp(1);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.Nodegroups{obj.myGroupIndex}.Nodes);
            
            obj.axisName = get_input(inputs, 'axis', "y");
            obj.side = get_input(inputs, 'side', "lower");
            obj.wallPosition = inputs.wallPosition;
            obj.pressurePenalty = get_input(inputs, 'pressurePenalty', 1e11);
            obj.activeTolerance = get_input(inputs, 'activeTolerance', 1e-12);
            obj.nodeArea = get_input(inputs, 'nodeArea', 1.0);
            obj.penalty = get_input(inputs, 'penalty', obj.pressurePenalty * obj.nodeArea);
            
            nnode = length(obj.mesh.Nodegroups{obj.myGroupIndex}.Nodes);
            obj.meanPressure = 0;
            obj.maxPressure = 0;
            obj.activeFraction = 0;
            obj.nodalPressure = zeros(nnode, 1);
            obj.nodalPenetration = zeros(nnode, 1);
        end
        
        function getKf(obj, physics, stp)
            if stp ~= obj.dx_Step
                return;
            end
            
            fprintf("        ContactPenalty (%s):", obj.myGroup)
            t = tic;
            
            nodes = obj.mesh.Nodegroups{obj.myGroupIndex}.Nodes;
            if obj.axisName == "x"
                comp = 1;
                dofType = obj.dofTypeIndices(1);
            else
                comp = 2;
                dofType = obj.dofTypeIndices(2);
            end
            dofs = obj.dofSpace.getDofIndices(dofType, nodes);
            q0 = obj.mesh.Nodes(nodes, comp);
            q = q0 + physics.StateVec{obj.dx_Step}(dofs);
            
            if obj.side == "lower"
                penetration = obj.wallPosition - q;
                signResidual = -1;
            else
                penetration = q - obj.wallPosition;
                signResidual = 1;
            end
            
            active = penetration >= -obj.activeTolerance;
            penetrationEff = max(penetration, 0);
            reaction = obj.penalty * penetrationEff;
            residual = signResidual * reaction;
            
            tangent = zeros(size(penetration));
            tangent(active) = obj.penalty;
            
            physics.fint{stp} = physics.fint{stp} + ...
                sparse(dofs, 0*dofs+1, residual, length(physics.fint{stp}), 1);
            physics.K{stp} = physics.K{stp} + ...
                sparse(dofs, dofs, tangent, length(physics.fint{stp}), length(physics.fint{stp}));
            
            pressure = obj.pressurePenalty * penetrationEff;
            obj.nodalPenetration = penetrationEff;
            obj.nodalPressure = pressure;
            obj.meanPressure = mean(pressure);
            obj.maxPressure = max(pressure);
            obj.activeFraction = nnz(active) / max(numel(active), 1);
            
            fprintf("  (%.2fs) meanPc=%.2fMPa maxPc=%.2fMPa active=%.2f\n", ...
                toc(t), obj.meanPressure/1e6, obj.maxPressure/1e6, obj.activeFraction);
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
