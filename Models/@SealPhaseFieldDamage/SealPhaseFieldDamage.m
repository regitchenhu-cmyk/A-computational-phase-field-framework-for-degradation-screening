classdef SealPhaseFieldDamage < BaseModel
    %SEALPHASEFIELDDAMAGE  Chemo-mechanically coupled degradation model
    %                      for within-study elastomer-seal screening.
    %
    % ===================== v11 — CLEAN REWRITE =====================
    %
    %  ─── 模型定位 ───
    %  相场 φ 代表「局部材料退化指标」（不是裂纹）:
    %    φ = 0: 材料完好
    %    φ = 1: 材料完全退化
    %
    %  失效判据：接触压力 < 液压压力 → 泄漏
    %    接触压力随 φ 下降: σ_contact ~ (1-φ)² · σ_initial
    %
    %  ─── 驱动力：多物理退化指标 ───
    %
    %  不再使用弹性应变能 Wplus 驱动相场！
    %  改用 CL(流体) + α(老化) + Df(疲劳) 的加权组合：
    %
    %    ψ_drive(x,t) = α_c · [ w_CL·CL + w_α·α + w_Df·Df ] · Gc/(2l)
    %
    %  这是物理上正确的：
    %    · 密封件退化由化学(CL,α)和疲劳(Df)驱动
    %    · 不由静态弹性能驱动（密封件设计成承受预压缩）
    %    · ψ_drive 的空间分布由 CL 扩散决定 → 左侧高、右侧低
    %    · 形成自然的退化梯度（从压力侧向非压力侧）
    %
    %  ─── 方程结构 ───
    %
    %  位移方程（含溶胀）:
    %    σ = g(φ)·D·(ε - α_sw·CL·I) + σ_residual
    %    g(φ) = (1-kmin)(1-φ)² + kmin
    %
    %  退化场方程（AT2 结构，退化驱动）:
    %    r = ∫ [(Gc/l + 2(1-kmin)H)·φ·N + Gc·l·∇φ·∇N
    %           - 2(1-kmin)H·N] dV
    %    H = max(H_old, ψ_drive)
    %
    %  耦合：
    %    CL → 溶胀 + Gc退化 + ψ_drive
    %    α  → Gc退化 + ψ_drive
    %    Df → Gc degradation + ψ_drive
    %    φ  → 刚度退化 → 接触压力下降 → 失效
    %    φ  → 增强扩散 D_eff (已在 FluidDiffusion 中)
    %
    % ═══════════════════════════════════════════════════════════

    properties
        mesh
        myName
        myGroup
        myGroupIndex
        dofSpace
        dofTypeIndices

        C10; D1; kmin; l; Gc
        constitutiveModel
        alpha_swell           % 溶胀系数
        agingDegrade          % 老化对 Gc 的退化系数
        fluidDegrade          % 流体对 Gc 的退化系数
        fatigueDegrade        % 疲劳对 Gc 的退化系数

        % 退化驱动权重
        alpha_coupling        % 总耦合强度
        degradeThreshold      % 退化驱动阈值
        w_CL                  % 流体权重
        w_aging               % 老化权重
        w_fatigue             % 疲劳权重
        w_pressure            % hydraulic pressure severity weight
        P_max                 % peak service pressure [Pa]
        P_ref                 % reference pressure [Pa]

        dx_Step; phi_step; CL_step; aging_step; fatigue_step

        doInit
        Hist; HistOld
        LFrac
        extentThreshold
        SpanningPath
        ContactPressureRatio  % 接触压力 / 液压压力
        T = 293.15;
    end

    methods
        function obj = SealPhaseFieldDamage(mesh, physics, inputs)
            obj.myName = "SealPhaseFieldDamage";
            disp("Initializing " + obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;

            [obj.dofTypeIndices, stp] = obj.dofSpace.getDofType({"dx","dy","phi","CL","alpha_a","D_f"});
            obj.dx_Step      = stp(1);
            obj.phi_step     = stp(3);
            obj.CL_step      = stp(4);
            obj.aging_step   = stp(5);
            obj.fatigue_step = stp(6);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));

            obj.C10 = inputs.C10;
            obj.D1  = inputs.D1;
            obj.kmin = inputs.kmin;
            if ~isscalar(obj.kmin) || obj.kmin < 0 || obj.kmin >= 1
                error('SealPhaseFieldDamage:InvalidResidualStiffness', ...
                    'kmin must be a scalar in the interval [0,1).');
            end
            obj.constitutiveModel = string(get_input(inputs, 'constitutiveModel', "linearized"));
            obj.l    = inputs.l;
            obj.Gc   = inputs.Gc;
            obj.alpha_swell  = inputs.alpha_swell;
            obj.agingDegrade = inputs.agingDegrade;
            obj.fluidDegrade = inputs.fluidDegrade;
            if isfield(inputs, 'fatigueDegrade')
                obj.fatigueDegrade = inputs.fatigueDegrade;
            else
                obj.fatigueDegrade = 0.0;
            end

            % 退化驱动参数
            obj.alpha_coupling = inputs.alpha_coupling;
            if isfield(inputs, 'degradeThreshold')
                obj.degradeThreshold = inputs.degradeThreshold;
            else
                obj.degradeThreshold = 0.0;
            end
            obj.w_CL           = inputs.w_CL;
            obj.w_aging        = inputs.w_aging;
            obj.w_fatigue      = inputs.w_fatigue;
            obj.w_pressure     = get_input(inputs, 'w_pressure', 0.0);
            obj.P_max          = get_input(inputs, 'P_max', 0.0);
            obj.P_ref          = get_input(inputs, 'P_ref', 35e6);

            nelem  = size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1);
            nipmax = obj.mesh.ipcount1D^2;
            obj.Hist    = zeros(nelem, nipmax);
            obj.HistOld = obj.Hist;
            obj.doInit  = true;
            obj.LFrac   = 0;
            obj.extentThreshold = get_input(inputs, 'extentThreshold', 0.35);
            obj.SpanningPath = false;
            obj.ContactPressureRatio = 1.0;
        end

        function OncePerStep(obj, physics, stp)
            if stp == obj.phi_step && obj.doInit
                % v11: 不预设种子裂纹, φ 从 0 自然演化
                % 但把初始 H 设为 0（退化驱动尚未累积）
                obj.doInit = false;
                fprintf("    [Init] φ starts from 0 (no seed crack)\n");
                fprintf("    [Init] Degradation-driven: α_c=%.1f, w_CL=%.2f, w_α=%.2f, w_Df=%.2f\n", ...
                    obj.alpha_coupling, obj.w_CL, obj.w_aging, obj.w_fatigue);
                fprintf("    [Init] H_crit = Gc/(2l) = %.2e J/m³\n", obj.Gc/(2*obj.l));
            end
        end

        function Commit(obj, ~, commit_type)
            if commit_type == "Pathdep"
                obj.HistOld = obj.Hist;
            end
        end

        function getKf(obj, physics, stp)
            mu     = 2*obj.C10;
            kappa  = 2/obj.D1;
            lambda = kappa - 2*mu/3;
            D_el = zeros(4,4);
            D_el(1,1) = lambda+2*mu; D_el(1,2) = lambda;      D_el(1,3) = lambda;
            D_el(2,1) = lambda;      D_el(2,2) = lambda+2*mu; D_el(2,3) = lambda;
            D_el(3,1) = lambda;      D_el(3,2) = lambda;      D_el(3,3) = lambda+2*mu;
            D_el(4,4) = mu;

            %% =====================================================
            %% Displacement step (含溶胀本征应变)
            %% =====================================================
            if stp == obj.dx_Step
                fprintf("        SealPFD (u):")
                t = tic;
                dofmatX=[]; dofmatY=[]; kmat=[]; fvec=[]; dofvec=[];
                SVec   = physics.StateVec;
                kmin_l = obj.kmin;
                asw    = obj.alpha_swell;

                for n_el = 1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                    [N, G, w]  = obj.mesh.getVals(obj.myGroupIndex, n_el);
                    dofsX  = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                    dofsY  = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                    dofsXY = [dofsX; dofsY];
                    dofsPhi = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
                    dofsCL  = obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), Elem_Nodes);
                    XY  = [SVec{obj.dx_Step}(dofsX); SVec{obj.dx_Step}(dofsY)];
                    PHI = SVec{obj.phi_step}(dofsPhi);
                    CL  = SVec{obj.CL_step}(dofsCL);

                    f_el = zeros(length(dofsXY), 1);
                    K_el = zeros(length(dofsXY));
                    for ip = 1:length(w)
                        B      = obj.getB(G(ip,:,:));
                        strain = B * XY;
                        ff     = min(max(N(ip,:)*PHI, 0), 1);
                        dam_fun = (1-kmin_l)*(1-ff)^2 + kmin_l;

                        % 溶胀本征应变
                        CL_ip = max(0, min(1, N(ip,:)*CL));
                        eps_sw = asw * CL_ip;
                        e_sw = [eps_sw; eps_sw; eps_sw; 0];

                        % 应力 = g(φ)·D·(ε - ε_sw)
                        strain_mech = strain - e_sw;
                        if obj.constitutiveModel == "neoHookean"
                            stress = dam_fun * obj.neoHookeanStress(G(ip,:,:), XY, eps_sw, mu, kappa);
                        else
                            stress = dam_fun * D_el * strain_mech;
                        end

                        f_el = f_el + B'*stress*w(ip);
                        K_el = K_el + B'*dam_fun*D_el*B*w(ip);
                    end
                    [dofmatxloc, dofmatyloc] = ndgrid(dofsXY, dofsXY);
                    dofmatX = [dofmatX; dofmatxloc(:)];
                    dofmatY = [dofmatY; dofmatyloc(:)];
                    kmat   = [kmat;   K_el(:)];
                    fvec   = [fvec;   f_el];
                    dofvec = [dofvec; dofsXY];
                end
                physics.fint{stp} = physics.fint{stp} + ...
                    sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{stp}), 1);
                physics.K{stp} = physics.K{stp} + ...
                    sparse(dofmatX, dofmatY, kmat, length(physics.fint{stp}), length(physics.fint{stp}));
                fprintf("  (%.2fs)\n", toc(t));
            end

            %% =====================================================
            %% Degradation field step (AT2 structure, degradation-driven)
            %%
            %%  ψ_drive = α_c · (w_CL·CL + w_α·α + w_Df·Df) · Gc/(2l)
            %%  H = max(H_old, ψ_drive)
            %%
            %%  r = ∫ [(Gc/l + 2(1-kmin)H)·φ·N + Gc·l·∇φ·∇N
            %%           - 2(1-kmin)H·N] dV
            %%  K = ∫ [(Gc/l + 2(1-kmin)H)·N·N + Gc·l·∇N·∇N] dV
            %%
            %%  Gc 本身也随退化降低:
            %%    Gc_eff = Gc · (1-f_CL·CL) · (1-f_α·α) · (1-f_Df·Df)
            %% =====================================================
            if stp == obj.phi_step
                fprintf("        SealPFD (phi):")
                t = tic;
                dofmatX=[]; dofmatY=[]; kmat=[]; fvec=[]; dofvec=[];
                SVec     = physics.StateVec;
                maxIP    = size(obj.Hist, 2);
                nelem    = size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1);
                HistOld_local = obj.HistOld;
                HNew_cell = cell(nelem, 1);
                Gc_g      = obj.Gc;
                l_g       = obj.l;
                fluidDg   = obj.fluidDegrade;
                agingDg   = obj.agingDegrade;
                fatigueDg = obj.fatigueDegrade;
                ac        = obj.alpha_coupling;
                dThresh   = obj.degradeThreshold;
                wCL       = obj.w_CL;
                wA        = obj.w_aging;
                wDf       = obj.w_fatigue;
                kmin_l    = obj.kmin;
                degradationSlope = 1 - kmin_l;
                pressureTerm = obj.w_pressure * min(1, max(0, obj.P_max / obj.P_ref));

                for n_el = 1:nelem
                    Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                    [N, G, w]  = obj.mesh.getVals(obj.myGroupIndex, n_el);
                    dofsPhi   = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
                    dofsCL    = obj.dofSpace.getDofIndices(obj.dofTypeIndices(4), Elem_Nodes);
                    dofsAging = obj.dofSpace.getDofIndices(obj.dofTypeIndices(5), Elem_Nodes);
                    dofsFat   = obj.dofSpace.getDofIndices(obj.dofTypeIndices(6), Elem_Nodes);

                    PHI     = SVec{obj.phi_step}(dofsPhi);
                    CL      = SVec{obj.CL_step}(dofsCL);
                    AGING   = SVec{obj.aging_step}(dofsAging);
                    FATIGUE = SVec{obj.fatigue_step}(dofsFat);

                    f_el    = zeros(length(dofsPhi), 1);
                    K_el    = zeros(length(dofsPhi));
                    H_local = zeros(1, maxIP);

                    for ip = 1:maxIP
                        Bgrad = squeeze(G(ip,:,:));

                        % 各退化变量在积分点的值
                        CL_ip  = max(0, min(1, N(ip,:)*CL));
                        ag_ip  = max(0, min(1, N(ip,:)*AGING));
                        Df_ip  = max(0, min(1, N(ip,:)*FATIGUE));

                        % Gc 退化
                        Gc_eff = Gc_g * (1 - fluidDg*CL_ip) * ...
                            (1 - agingDg*ag_ip) * (1 - fatigueDg*Df_ip);
                        Gc_eff = max(Gc_eff, 0.02*Gc_g);

                        % 退化驱动力（多物理加权）
                        degrade_index = wCL*CL_ip + wA*ag_ip + wDf*Df_ip + pressureTerm;
                        if dThresh > 0
                            degrade_index = max(0, degrade_index - dThresh) / ...
                                max(1e-12, 1 - dThresh);
                        end
                        psi_drive = ac * degrade_index * Gc_eff / (2*l_g);

                        % Miehe 单调历史
                        H_new       = max(HistOld_local(n_el, ip), psi_drive);
                        H_local(ip) = H_new;

                        % AT2 方程（标准结构）
                        coef_lin = Gc_eff/l_g + 2*degradationSlope*H_new;
                        coef_grd = Gc_eff * l_g;

                        f_el = f_el + w(ip)*( ...
                            coef_lin*(N(ip,:)'*N(ip,:))*PHI + ...
                            coef_grd*(Bgrad*Bgrad')*PHI - ...
                            2*degradationSlope*H_new*N(ip,:)' );

                        K_el = K_el + w(ip)*( ...
                            coef_lin*(N(ip,:)'*N(ip,:)) + ...
                            coef_grd*(Bgrad*Bgrad') );
                    end

                    HNew_cell{n_el} = H_local;
                    [dofmatxloc, dofmatyloc] = ndgrid(dofsPhi, dofsPhi);
                    dofmatX = [dofmatX; dofmatxloc(:)];
                    dofmatY = [dofmatY; dofmatyloc(:)];
                    kmat   = [kmat;   K_el(:)];
                    fvec   = [fvec;   f_el];
                    dofvec = [dofvec; dofsPhi];
                end

                physics.fint{stp} = physics.fint{stp} + ...
                    sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint{stp}), 1);
                physics.K{stp} = physics.K{stp} + ...
                    sparse(dofmatX, dofmatY, kmat, length(physics.fint{stp}), length(physics.fint{stp}));

                HNew = zeros(size(obj.Hist));
                for n_el = 1:nelem
                    if ~isempty(HNew_cell{n_el})
                        HNew(n_el,:) = HNew_cell{n_el};
                    end
                end
                obj.Hist = HNew;

                % φ 统计和接触压力比
                allNodes = obj.mesh.GetAllNodesForGroup(obj.myGroupIndex);
                PhiDofs  = obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), allNodes);
                phi_nodes = physics.StateVec{obj.phi_step}(PhiDofs);
                phi_max   = max(phi_nodes);
                phi_avg   = mean(phi_nodes);

                % 接触压力比: 用右侧面（x=max）节点的 φ 估算
                x_all = obj.mesh.Nodes(allNodes, 1);
                x_max = max(x_all);
                right_mask = (x_all > x_max - 1e-6);
                if any(right_mask)
                    phi_right = mean(phi_nodes(right_mask));
                    obj.ContactPressureRatio = (1-phi_right)^2;
                else
                    obj.ContactPressureRatio = (1-phi_avg)^2;
                end

                % Left-connected thresholded degradation extent. This is
                % deliberately not called a crack length. A breadth-first
                % search prevents isolated or uniformly sub-threshold nodes
                % from being reported as a macroscopic path.
                y_all = obj.mesh.Nodes(allNodes, 2);
                [x_vals, ~, ix] = unique(round(x_all, 12));
                [y_vals, ~, iy] = unique(round(y_all, 12));
                open_grid = false(numel(y_vals), numel(x_vals));
                for nn = 1:numel(phi_nodes)
                    open_grid(iy(nn), ix(nn)) = open_grid(iy(nn), ix(nn)) || ...
                        phi_nodes(nn) >= obj.extentThreshold;
                end
                visited = false(size(open_grid));
                queue = zeros(numel(open_grid), 2);
                head = 1; tail = 0;
                for jj = 1:size(open_grid, 1)
                    if open_grid(jj, 1)
                        tail = tail + 1;
                        queue(tail,:) = [jj, 1];
                        visited(jj, 1) = true;
                    end
                end
                while head <= tail
                    jj = queue(head, 1);
                    ii = queue(head, 2);
                    head = head + 1;
                    neigh = [jj-1 ii; jj+1 ii; jj ii-1; jj ii+1];
                    for kk = 1:4
                        j2 = neigh(kk, 1); i2 = neigh(kk, 2);
                        if j2 >= 1 && j2 <= size(open_grid,1) && ...
                                i2 >= 1 && i2 <= size(open_grid,2) && ...
                                open_grid(j2,i2) && ~visited(j2,i2)
                            tail = tail + 1;
                            queue(tail,:) = [j2, i2];
                            visited(j2,i2) = true;
                        end
                    end
                end
                [~, connected_cols] = find(visited);
                if isempty(connected_cols)
                    obj.LFrac = 0;
                else
                    obj.LFrac = max(x_vals(connected_cols)) - min(x_vals);
                end
                obj.SpanningPath = any(visited(:, end));

                fprintf("  (%.2fs) maxH=%.2e  φ_avg=%.3f  φ_max=%.3f  CP_ratio=%.3f  connected_extent=%.2fmm  span=%d\n", ...
                    toc(t), max(HNew(:)), phi_avg, phi_max, obj.ContactPressureRatio, ...
                    obj.LFrac*1e3, obj.SpanningPath);
            end
        end

        function B = getB(~, grads)
            cp_count = size(grads, 2);
            B = zeros(4, cp_count*2);
            for ii = 1:cp_count
                B(1, ii)          = grads(1,ii,1);
                B(4, ii)          = grads(1,ii,2);
                B(2, ii+cp_count) = grads(1,ii,2);
                B(4, ii+cp_count) = grads(1,ii,1);
            end
        end

        function stress = neoHookeanStress(~, grads, XY, eps_sw, mu, kappa)
            cp_count = size(grads, 2);
            ux = XY(1:cp_count);
            uy = XY(cp_count+1:end);
            dNdx = squeeze(grads(1,:,1))';
            dNdy = squeeze(grads(1,:,2))';

            F2 = [1 + ux' * dNdx, ux' * dNdy; ...
                  uy' * dNdx,     1 + uy' * dNdy];
            swell = max(1e-6, 1 + eps_sw);
            F2 = F2 / swell;

            F = eye(3);
            F(1:2,1:2) = F2;
            J = max(det(F), 1e-8);
            b = F * F';

            sigma = (mu / J) * (b - eye(3)) + kappa * (J - 1) * eye(3);
            stress = [sigma(1,1); sigma(2,2); sigma(3,3); sigma(1,2)];
        end

        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
            hasInfo = false; provided = [];
            if (var=="stresses"||var=="sxx"||var=="syy"||var=="szz"||var=="sxy"||var=="sh")
                hasInfo = true;
                mu=2*obj.C10; kappa=2/obj.D1; lambda=kappa-2*mu/3;
                D_el=zeros(4); D_el(1,1)=lambda+2*mu; D_el(1,2)=lambda; D_el(1,3)=lambda;
                D_el(2,1)=lambda; D_el(2,2)=lambda+2*mu; D_el(2,3)=lambda;
                D_el(3,1)=lambda; D_el(3,2)=lambda; D_el(3,3)=lambda+2*mu; D_el(4,4)=mu;
                if var=="stresses"
                    provided=zeros(length(elems),obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount,4);
                else
                    provided=zeros(length(elems),obj.mesh.Elementgroups{obj.myGroupIndex}.ShapeFunc.ipcount);
                end
                for el = 1:length(elems)
                    Elem_Nodes=obj.mesh.getNodes(obj.myGroupIndex,elems(el));
                    [N,G,w]=obj.mesh.getVals(obj.myGroupIndex,elems(el));
                    dofsX=obj.dofSpace.getDofIndices(obj.dofTypeIndices(1),Elem_Nodes);
                    dofsY=obj.dofSpace.getDofIndices(obj.dofTypeIndices(2),Elem_Nodes);
                    dofsPhi=obj.dofSpace.getDofIndices(obj.dofTypeIndices(3),Elem_Nodes);
                    dofsCL=obj.dofSpace.getDofIndices(obj.dofTypeIndices(4),Elem_Nodes);
                    XY=[physics.StateVec{obj.dx_Step}(dofsX);physics.StateVec{obj.dx_Step}(dofsY)];
                    PHI=physics.StateVec{obj.phi_step}(dofsPhi);
                    CL=physics.StateVec{obj.CL_step}(dofsCL);
                    for ip=1:length(w)
                        ff=min(max(N(ip,:)*PHI,0),1);
                        dam=(1-obj.kmin)*(1-ff)^2+obj.kmin;
                        B=obj.getB(G(ip,:,:));
                        str=B*XY;
                        CL_ip=max(0,min(1,N(ip,:)*CL));
                        esw=obj.alpha_swell*CL_ip;
                        str_m=str-[esw;esw;esw;0];
                        stress=dam*D_el*str_m;
                        if loc=="Interior"
                            switch var
                                case "stresses"; provided(el,ip,:)=stress;
                                case "sxx"; provided(el,ip)=stress(1);
                                case "syy"; provided(el,ip)=stress(2);
                                case "szz"; provided(el,ip)=stress(3);
                                case "sxy"; provided(el,ip)=stress(4);
                                case "sh"; provided(el,ip)=(stress(1)+stress(2)+stress(3))/3;
                            end
                        end
                    end
                end
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
