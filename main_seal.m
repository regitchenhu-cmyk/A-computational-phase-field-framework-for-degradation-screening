function [physics, tvec, results] = main_seal(scenario, overrides)
    %MAIN_SEAL Coupled degradation-screening simulation for a hydraulic seal.
    %
    % Simulates prescribed comparative operating cases by coupling:
    %   1. AT2-regularized diffuse degradation
    %   2. Fatigue damage (cyclic pressure loading)
    %   3. Hydraulic fluid diffusion (swelling & property degradation)
    %   4. Effective aging-state evolution
    %
    % Usage:
    %   [physics, tvec, results] = main_seal(1)  % Standard O-ring scenario
    %   [physics, tvec, results] = main_seal(2)  % High-temperature scenario
    %   [physics, tvec, results] = main_seal(3)  % High-pressure scenario
    %
    % Scenario definitions:
    %   1: Standard service - 21 MPa, 80°C, normal cycling
    %   2: High-temperature - 21 MPa, 120°C, accelerated aging
    %   3: High-pressure    - 35 MPa, 80°C, aggressive fatigue
    %
    % IMPORTANT: the direct-call defaults below are retained for historical
    % demonstrations and are not the canonical paper parameter set. Recreate
    % paper cases through run_scientifically_corrected_case_20260823, which
    % applies the prescribed values in calibrated_degradation_options.m.
    
    if nargin < 1
        scenario = 1;
    end
    if nargin < 2 || isempty(overrides)
        overrides = struct();
    end
    
    %% Scenario configuration
    switch scenario
        case 1  % Standard service conditions
            P_max = 21e6;       % 21 MPa (3000 psi typical aircraft hydraulic)
            T_service = 353.15; % 80°C
            f_cycle = 0.05;     % 0.05 Hz (slow pressure cycling)
            sname = "Standard_21MPa_80C";
        case 2  % High-temperature conditions
            P_max = 21e6;
            T_service = 393.15; % 120°C (near engine/hot zone)
            f_cycle = 0.05;
            sname = "HighTemp_21MPa_120C";
        case 3  % High-pressure conditions
            P_max = 35e6;       % 35 MPa high-pressure aircraft hydraulic case
            % === CYCLIC PRESSURE (uncomment to enable) ===
            % freq = 1; P_min = 0;
            % p_h = P_min + (P_max - P_min) * (0.5 + 0.5*sin(2*pi*freq*t));
            T_service = 353.15;
            f_cycle = 0.1;      % Faster cycling
            sname = "HighPress_35MPa_80C";
    end
    if isfield(overrides, 'P_max'), P_max = overrides.P_max; end
    if isfield(overrides, 'T_service'), T_service = overrides.T_service; end
    if isfield(overrides, 'f_cycle'), f_cycle = overrides.f_cycle; end
    if isfield(overrides, 'sname'), sname = string(overrides.sname); end
    if isfield(overrides, 'nameSuffix'), sname = sname + "_" + string(overrides.nameSuffix); end
    
    maxThreads = get_option(overrides, 'maxThreads', 16);
    maxNumCompThreads(maxThreads);
    
    saveRoot = string(get_option(overrides, 'saveRoot', "./Results_Seal"));
    savefolder = saveRoot + "/" + sname;
    mkdir(savefolder);
    savefolder = savefolder + "/";
    
    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))
    
    useParpool = get_option(overrides, 'useParpool', true);
    if useParpool && isempty(gcp('nocreate'))
        parpool('threads')
    end
    
    tmr = tic;
    
    %% Check for restart
    % By default, always start fresh. Set force_restart = true to 
    % continue from a previous run (requires compatible class definitions).
    force_restart = get_option(overrides, 'forceRestart', false);
    
    Files = dir(fullfile(savefolder, '*.mat'));
    nmax_existing = length(Files);
    
    if force_restart && nmax_existing > 0
        % Find the highest numeric step file
        restart_num = 0;
        for fi = 1:length(Files)
            [~, fname, ~] = fileparts(Files(fi).name);
            num = str2double(fname);
            if ~isnan(num) && num > restart_num
                restart_num = num;
            end
        end
        if restart_num > 0
            restart = true;
            fprintf('Restarting from step %d\n', restart_num);
        else
            restart = false;
        end
    else
        restart = false;
        % Clean out any old results to prevent class conflicts
        cleanResults = get_option(overrides, 'cleanResults', true);
        if cleanResults && nmax_existing > 0
            fprintf('Clearing %d old result files in %s\n', nmax_existing, savefolder);
            delete(fullfile(savefolder, '*.mat'));
        end
    end
    
    if restart == false
        %% =============================================
        %%  MESH DEFINITION - Seal Cross-Section
        %% =============================================
        % Rectangular cross-section of an O-ring seal
        % (simplification of actual circular profile)
        mesh_in.type    = "Seal";
        mesh_in.Nx      = 40;        % Elements in radial direction
        mesh_in.Ny      = 56;        % Elements in axial direction
        mesh_in.Lx      = 3.5e-3;    % Seal width = 3.5 mm (radial)
        mesh_in.Ly      = 5.0e-3;    % Seal height = 5.0 mm (axial)
        mesh_in.ipcount1D = 3;        % 3x3 Gauss points
        mesh_in.zeroWeight = false;
        mesh_in.grooveDepth = 1.0e-3; % 1mm groove constraint
        if isfield(overrides, 'mesh')
            mesh_in = merge_struct(mesh_in, overrides.mesh);
        end
        % Element size: hx = 3.5/40 = 0.0875 mm, hy = 5/56 = 0.0893 mm
        % Phase-field length scale chosen as l = 2h ~ 0.18 mm
        
        %% =============================================
        %%  DOF DEFINITION
        %% =============================================
        %  dx, dy  : Displacement (mechanical equilibrium)
        %  phi     : Phase-field damage variable
        %  CL      : Hydraulic fluid concentration
        %  alpha_a : Aging degradation variable
        %  D_f     : Fatigue damage variable
        dofs_in.dofs = {"dx","dy","phi","CL","alpha_a","D_f"};
        dofs_in.Step = [  2,   2,   1,   3,     4,      4  ];
        %  Step 1: Phase-field evolution
        %  Step 2: Mechanical equilibrium
        %  Step 3: Fluid diffusion
        %  Step 4: Aging + Fatigue (coupled in same step)
        
        %% =============================================
        %%  PHYSICS MODELS
        %% =============================================
        % Historical demonstration parameters were numerically tuned such that:
        %  - Mechanical strain energy at crack tip is ~10*W0 → fatigue active
        %  - Pressure-driven Wplus at tip > Gc/(2l) → phase-field grows
        %  - Aging timescale 1/k_eff ≈ 30 days at 80°C
        %  - Fluid diffusion timescale L²/D ≈ 30 days for L=3.5mm
        
        %% Model 1: Chemo-mechanically coupled degradation (v11)
        physics_in{1}.type           = "SealPhaseFieldDamage";
        physics_in{1}.Egroup         = "Internal";
        physics_in{1}.C10            = 1.5e6;    % Mooney-Rivlin C10 [Pa]
        physics_in{1}.constitutiveModel = string(get_option(overrides, 'constitutiveModel', "linearized"));
        physics_in{1}.D1             = 2e-9;     % 近不可压缩
        physics_in{1}.kmin           = 5e-3;     % 残余刚度
        physics_in{1}.l              = 0.30e-3;  % 正则化长度 [m]
        physics_in{1}.Gc             = 200;      % 断裂韧性 [J/m²]
        physics_in{1}.alpha_swell    = 0.04;     % 溶胀系数
        physics_in{1}.agingDegrade   = 0.7;      % α 退化 Gc 的系数
        physics_in{1}.fluidDegrade   = 0.5;      % CL 退化 Gc 的系数
        physics_in{1}.fatigueDegrade = 0.35;     % D_f 退化 Gc 的系数
        % 退化驱动参数
        %   ψ_drive = α_c · (w_CL·CL + w_α·α + w_Df·Df) · Gc_eff/(2l)
        %   当加权退化指标 > 1/α_c 时 φ 开始显著增长
        physics_in{1}.alpha_coupling = 2.0;      % 耦合强度
        physics_in{1}.w_CL           = 0.4;      % 流体权重
        physics_in{1}.w_aging        = 0.4;      % 老化权重
        physics_in{1}.w_fatigue      = 0.2;      % 疲劳权重
        physics_in{1}.w_pressure     = 0.0;
        physics_in{1}.P_max          = P_max;
        physics_in{1}.P_ref          = 35e6;
        if isfield(overrides, 'alpha_coupling'), physics_in{1}.alpha_coupling = overrides.alpha_coupling; end
        if isfield(overrides, 'degradeThreshold'), physics_in{1}.degradeThreshold = overrides.degradeThreshold; end
        if isfield(overrides, 'w_CL'), physics_in{1}.w_CL = overrides.w_CL; end
        if isfield(overrides, 'w_aging'), physics_in{1}.w_aging = overrides.w_aging; end
        if isfield(overrides, 'w_fatigue'), physics_in{1}.w_fatigue = overrides.w_fatigue; end
        if isfield(overrides, 'w_pressure'), physics_in{1}.w_pressure = overrides.w_pressure; end
        if isfield(overrides, 'P_ref'), physics_in{1}.P_ref = overrides.P_ref; end
        if isfield(overrides, 'fluidDegrade'), physics_in{1}.fluidDegrade = overrides.fluidDegrade; end
        if isfield(overrides, 'agingDegrade'), physics_in{1}.agingDegrade = overrides.agingDegrade; end
        if isfield(overrides, 'fatigueDegrade'), physics_in{1}.fatigueDegrade = overrides.fatigueDegrade; end
        if isfield(overrides, 'l'), physics_in{1}.l = overrides.l; end
        if isfield(overrides, 'Gc'), physics_in{1}.Gc = overrides.Gc; end
        if isfield(overrides, 'C10'), physics_in{1}.C10 = overrides.C10; end
        if isfield(overrides, 'D1'), physics_in{1}.D1 = overrides.D1; end
        physics_in{1}.extentThreshold = get_option(overrides, 'degradationExtentThreshold', 0.35);
        
        %% Model 2: Hydraulic fluid diffusion (simple Fick's law)
        physics_in{2}.type      = "FluidDiffusion";
        physics_in{2}.Egroup    = "Internal";
        % Calibrated so diffusion length √(Dt) ~ Lx/2 in 60 days
        %   D ≈ (Lx/2)^2 / 60d = (1.75e-3)^2/(5.2e6) = 6e-13
        physics_in{2}.DL        = 6e-13;
        physics_in{2}.kmin      = 5e-3;
        physics_in{2}.DL_crack  = 1e-9;
        if isfield(overrides, 'DL'), physics_in{2}.DL = overrides.DL; end
        
        squeeze_ratio = 0.08;  % PATCHED: 8% pre-compression to create tensile zones: pressure alone drives the seal
                               %   (squeeze added 40+ MPa of pre-compression
                               %    that overwhelmed the 21 MPa hydraulic
                               %    pressure-driven opening force)
        if isfield(overrides, 'squeeze_ratio'), squeeze_ratio = overrides.squeeze_ratio; end
        
        useContactPenalty = get_option(overrides, 'useContactPenalty', true);
        contactPressure0 = get_option(overrides, 'contactPressure0', 1.2*P_max);
        contactPenaltyPressure = get_option(overrides, 'contactPenaltyPressure', ...
            contactPressure0 / max(squeeze_ratio * mesh_in.Ly, eps));
        contactTol = get_option(overrides, 'contactActiveTolerance', 1e-10);
        topNodeArea = get_option(overrides, 'topContactNodeArea', mesh_in.Lx/(2*mesh_in.Nx));
        bottomNodeArea = get_option(overrides, 'bottomContactNodeArea', mesh_in.Lx/(2*mesh_in.Nx));
        rightNodeArea = get_option(overrides, 'rightContactNodeArea', mesh_in.Ly/(2*mesh_in.Ny));
        
        if useContactPenalty
            %% Model 3: Bottom boundary - groove floor contact
            physics_in{3}.type   = "ContactPenalty";
            physics_in{3}.Ngroup = "Bottom";
            physics_in{3}.axis   = "y";
            physics_in{3}.side   = "lower";
            physics_in{3}.wallPosition = 0;
            physics_in{3}.pressurePenalty = contactPenaltyPressure;
            physics_in{3}.nodeArea = bottomNodeArea;
            physics_in{3}.activeTolerance = contactTol;
            
            %% Model 4: Right boundary - groove wall contact
            physics_in{4}.type   = "ContactPenalty";
            physics_in{4}.Ngroup = "Right";
            physics_in{4}.axis   = "x";
            physics_in{4}.side   = "upper";
            physics_in{4}.wallPosition = mesh_in.Lx;
            physics_in{4}.pressurePenalty = contactPenaltyPressure;
            physics_in{4}.nodeArea = rightNodeArea;
            physics_in{4}.activeTolerance = contactTol;
            
            %% Model 5: Top boundary - mating surface contact
            physics_in{5}.type   = "ContactPenalty";
            physics_in{5}.Ngroup = "Top";
            physics_in{5}.axis   = "y";
            physics_in{5}.side   = "upper";
            physics_in{5}.wallPosition = mesh_in.Ly - squeeze_ratio * mesh_in.Ly;
            physics_in{5}.pressurePenalty = contactPenaltyPressure;
            physics_in{5}.nodeArea = topNodeArea;
            physics_in{5}.activeTolerance = contactTol;
        else
            %% Model 3: Bottom boundary - groove floor (dy = 0)
            physics_in{3}.type   = "Constrainer";
            physics_in{3}.Ngroup = "Bottom";
            physics_in{3}.dofs   = {"dy"};
            physics_in{3}.conVal = [0];
            
            %% Model 4: Right boundary - groove wall (dx = 0)
            physics_in{4}.type   = "Constrainer";
            physics_in{4}.Ngroup = "Right";
            physics_in{4}.dofs   = {"dx"};
            physics_in{4}.conVal = [0];
            
            %% Model 5: Top boundary - mating surface (dy = -squeeze)
            physics_in{5}.type   = "Constrainer";
            physics_in{5}.Ngroup = "Top";
            physics_in{5}.dofs   = {"dy"};
            physics_in{5}.conVal = [-squeeze_ratio * mesh_in.Ly];
        end
        
        %% Model 6: Pressure load on high-pressure side (left)
        physics_in{6}.type       = "PressureLoad";
        physics_in{6}.Egroup     = "Left";
        physics_in{6}.pressure   = P_max;
        physics_in{6}.direction  = "normal";
        physics_in{6}.P_profile  = "ramp";
        physics_in{6}.P_min      = 0.5e6;
        physics_in{6}.f_cycle    = f_cycle;
        physics_in{6}.t_ramp     = 1000;
        if isfield(overrides, 'P_min'), physics_in{6}.P_min = overrides.P_min; end
        
        %% Model 7: Aging degradation
        % Calibrated so alpha=~0.3 (mild aging) at end of 60-day sim:
        %   k_eff = -ln(0.7)/(60 days) ~ 7e-8 /s
        %   k_eff = k_age*33 at 80C  =>  k_age ~ 2e-9
        physics_in{7}.type      = "AgingDegradation";
        physics_in{7}.Egroup    = "Internal";
        physics_in{7}.k_age     = 2e-9;
        physics_in{7}.n_age     = 1.0;
        physics_in{7}.Ea_age    = 50e3;
        physics_in{7}.D_age     = 5e-12;        % slower diffusion of aging agent
        physics_in{7}.phiDiffusionFactor = get_option(overrides, 'agingPhiDiffusionFactor', 0.0);
        physics_in{7}.T_service = T_service;
        if isfield(overrides, 'k_age'), physics_in{7}.k_age = overrides.k_age; end
        
        %% Model 8: Fatigue damage (pointwise update + L2 projection)
        % W0 calibrated so that with tip Wplus ~ 1.5 MPa, Df reaches
        % ~0.5 over 30 days (130k cycles at f=0.05 Hz).
        physics_in{8}.type        = "FatigueDamage";
        physics_in{8}.Egroup      = "Internal";
        physics_in{8}.W0          = 1.5e8;      % J/m³  (effective fatigue threshold)
        physics_in{8}.beta        = 2.0;
        physics_in{8}.P_max       = P_max;
        physics_in{8}.P_min       = 0.5e6;
        physics_in{8}.P_ref       = 21e6;
        physics_in{8}.pressureExponent = 0.0;
        physics_in{8}.energyRateScale = 1.0;
        physics_in{8}.contactFatigueCoeff = 0.25;
        physics_in{8}.contactFatigueExponent = 1.0;
        physics_in{8}.contactPressure0 = 1.2 * P_max;
        physics_in{8}.contactPressureRef = P_max;
        physics_in{8}.contactDecayLength = 0.5e-3;
        physics_in{8}.baseRate    = 0.0;
        physics_in{8}.maxRate     = 1e-5;
        physics_in{8}.f_cycle     = f_cycle;
        physics_in{8}.kmin        = 5e-3;
        if isfield(overrides, 'W0'), physics_in{8}.W0 = overrides.W0; end
        if isfield(overrides, 'beta'), physics_in{8}.beta = overrides.beta; end
        if isfield(overrides, 'fatiguePressureExponent'), physics_in{8}.pressureExponent = overrides.fatiguePressureExponent; end
        if isfield(overrides, 'fatigueEnergyRate'), physics_in{8}.energyRateScale = overrides.fatigueEnergyRate; end
        if isfield(overrides, 'contactFatigueCoeff'), physics_in{8}.contactFatigueCoeff = overrides.contactFatigueCoeff; end
        if isfield(overrides, 'contactFatigueExponent'), physics_in{8}.contactFatigueExponent = overrides.contactFatigueExponent; end
        if isfield(overrides, 'contactPressure0'), physics_in{8}.contactPressure0 = overrides.contactPressure0; end
        if isfield(overrides, 'contactPressureRef'), physics_in{8}.contactPressureRef = overrides.contactPressureRef; end
        if isfield(overrides, 'contactDecayLength'), physics_in{8}.contactDecayLength = overrides.contactDecayLength; end
        if isfield(overrides, 'fatigueP_ref'), physics_in{8}.P_ref = overrides.fatigueP_ref; end
        if isfield(overrides, 'fatigueBaseRate'), physics_in{8}.baseRate = overrides.fatigueBaseRate; end
        if isfield(overrides, 'fatigueMaxRate'), physics_in{8}.maxRate = overrides.fatigueMaxRate; end
        if isfield(overrides, 'P_min'), physics_in{8}.P_min = overrides.P_min; end
        
        %% Model 9: Left-bottom corner fix (dx = 0) 
        physics_in{9}.type   = "Constrainer";
        physics_in{9}.Ngroup = "LeftBottom";
        physics_in{9}.dofs   = {"dx"};
        physics_in{9}.conVal = [0];
        
        %% Model 10: Fluid BC - left boundary exposed to hydraulic fluid
        physics_in{10}.type   = "Constrainer";
        physics_in{10}.Ngroup = "Left";
        physics_in{10}.dofs   = {"CL"};
        physics_in{10}.conVal = [1.0];
        if isfield(overrides, 'CL_boundary'), physics_in{10}.conVal = overrides.CL_boundary; end
        
        %% Model 11: Aging BC - prescribed right-surface aging extent
        % Lower the surface aging level so the field has room to evolve;
        % the bulk reaction k_age*(1-alpha) generates aging from inside too.
        physics_in{11}.type   = "Constrainer";
        physics_in{11}.Ngroup = "Right";
        physics_in{11}.dofs   = {"alpha_a"};
        physics_in{11}.conVal = [0.3];   % prescribed surface aging extent
        if isfield(overrides, 'alpha_boundary'), physics_in{11}.conVal = overrides.alpha_boundary; end
        
        %% Model 12: Leakage estimator - Darcy flow through damaged/contact-open paths
        Q_allow = get_option(overrides, 'Q_allow', 1e-12);
        physics_in{12}.type      = "LeakagePostProcessor";
        physics_in{12}.Egroup    = "Internal";
        physics_in{12}.P_high    = P_max;
        physics_in{12}.P_low     = get_option(overrides, 'P_low', 0.101e6);
        physics_in{12}.mu        = get_option(overrides, 'fluidViscosity', 0.046);
        physics_in{12}.k0        = get_option(overrides, 'leakage_k0', 1e-22);
        physics_in{12}.k_phi     = get_option(overrides, 'leakage_k_phi', 1e-14);
        physics_in{12}.phiExponent = get_option(overrides, 'leakagePhiExponent', 4.0);
        physics_in{12}.aperture0 = get_option(overrides, 'leakageAperture0', 1e-9);
        physics_in{12}.aperturePhi = get_option(overrides, 'leakageAperturePhi', 5e-6);
        physics_in{12}.apertureExponent = get_option(overrides, 'leakageApertureExponent', 2.0);
        physics_in{12}.contactChi = get_option(overrides, 'leakageContactChi', 3.0);
        physics_in{12}.contactDecayLength = get_option(overrides, 'contactDecayLength', 0.5e-3);
        physics_in{12}.contactPressure0 = get_option(overrides, 'contactPressure0', 1.2*P_max);
        physics_in{12}.contactPressureRef = get_option(overrides, 'contactPressureRef', P_max);
        physics_in{12}.outOfPlaneWidth = get_option(overrides, 'outOfPlaneWidth', 1.0);
        physics_in{12}.connectedPhiThreshold = get_option(overrides, 'connectedPhiThreshold', 0.35);
        physics_in{12}.connectedLeakinessThreshold = get_option(overrides, 'connectedLeakinessThreshold', 1e-3);
        
        %% =============================================
        %%  SOLVER SETTINGS
        %% =============================================
        solver_in.maxIt = 80;
        solver_in.Conv = 1e-5;
        solver_in.tiny = 1e-3;
        solver_in.linesearch = false;
        solver_in.linesearchLims = [0.1 1];
        solver_in.OuterLoops = 8;
        
        %% =============================================
        %%  INITIALIZATION
        %% =============================================
        mesh = Mesh(mesh_in);
        mesh.check();
        
        physics = Physics(mesh, physics_in, dofs_in);
        
        dt = get_option(overrides, 'initialDt', 30);  % Initial time step [s]
        physics.time = 0;
        
        solver = Solver(physics, solver_in);
        tvec = 0;
        
        % Result tracking arrays
        results.CL_avg   = 0;
        results.CL_max   = 0;
        results.LFrac    = 0;
        results.degradation_connected = false;
        results.phi_avg = 0;
        results.phi_p99 = 0;
        results.alpha_p99 = 0;
        results.Df_p99 = 0;
        results.screening_index = 0;
        results.alpha_avg = 0;
        results.alpha_max = 0;
        results.Df_max   = 0;
        results.totalCycles = 0;
        results.leakage_Q = 0;
        results.leakage_Q_poiseuille = 0;
        results.leakage_kmax = physics.models{12}.maxPermeability;
        results.leakage_aperture_max = physics.models{12}.maxAperture;
        results.leakage_path_cost = physics.models{12}.pathCost;
        results.leakage_path_mean_aperture = physics.models{12}.pathMeanAperture;
        results.leakage_path_min_aperture = physics.models{12}.pathMinAperture;
        results.leakage_connected = false;
        results.Q_allow = Q_allow;
        results.leakage_failed_darcy = false;
        results.leakage_failed_poiseuille = false;
        results.leakage_time_darcy = NaN;
        results.leakage_time_poiseuille = NaN;
        results.contactPressureMean = physics.models{12}.contactPressureMean;
        results.contactClosureMean = physics.models{12}.contactClosureMean;
        contactMetrics = collect_contact_metrics(physics);
        results.contact_top_mean = contactMetrics.topMean;
        results.contact_bottom_mean = contactMetrics.bottomMean;
        results.contact_right_mean = contactMetrics.rightMean;
        results.contact_top_max = contactMetrics.topMax;
        results.contact_bottom_max = contactMetrics.bottomMax;
        results.contact_right_max = contactMetrics.rightMax;
        
        n_max = get_option(overrides, 'nMax', 800);              % maximum time steps
        tmax  = get_option(overrides, 'tmaxDays', 90) * 86400;   % maximum time [s]
        
        startstep = 1;
    else
        filename = savefolder + string(restart_num);
        load(filename, "mesh","physics","solver","dt","tvec","results","n_max","tmax");
        n_max = get_option(overrides, 'nMax', n_max);
        tmax = get_option(overrides, 'tmaxDays', tmax / 86400) * 86400;
        startstep = restart_num + 1;
    end
    maxDt = get_option(overrides, 'maxDt', 6*3600);
    plotEvery = get_option(overrides, 'plotEvery', inf);
    saveEvery = get_option(overrides, 'saveEvery', 50);
    verboseEvery = get_option(overrides, 'verboseEvery', 30);
    leakageModelIndex = find_model(physics, "LeakagePostProcessor");
    if ~isfield(results, 'Q_allow'), results.Q_allow = get_option(overrides, 'Q_allow', 1e-12); end
    if ~isfield(results, 'leakage_failed_darcy'), results.leakage_failed_darcy = false(size(tvec)); end
    if ~isfield(results, 'leakage_failed_poiseuille'), results.leakage_failed_poiseuille = false(size(tvec)); end
    if ~isfield(results, 'leakage_time_darcy'), results.leakage_time_darcy = NaN; end
    if ~isfield(results, 'leakage_time_poiseuille'), results.leakage_time_poiseuille = NaN; end
    
    %% =============================================
    %%  TIME STEPPING LOOP
    %% =============================================
    for tstep = startstep:n_max
        disp("===========================================");
        disp("Step: " + string(tstep));
        disp("Time: " + string(physics.time) + " s (" + ...
            string(physics.time/3600) + " hours, " + ...
            string(physics.time/3600/24) + " days)");
        
        % Adaptive time stepping
        % - Steps 1-50: small dt ~ 30s..2min, resolve mechanical/initial
        % - Steps 50-200: dt ~ minutes..hours, resolve diffusion onset
        % - Steps >200: dt ~ hours..6h, resolve aging+fatigue
        if tstep <= 30
            physics.dt = dt;                      % first 30 steps fixed: 30 s
        else
            physics.dt = dt * 1.04^(min(tstep-30, 200));
        end
        physics.dt = min(physics.dt, maxDt);      % Max time step
        disp("dTime: " + string(physics.dt) + " s (" + ...
            string(physics.dt/3600) + " hours)");
        
        % Pass current hydraulic pressure to phase-field model so it
        % can apply pressure body force inside the diffuse crack zone.
        % Use values stored on the PressureLoad model (model 6).
        t_now = physics.time + physics.dt;
        t_ramp = physics.models{6}.t_ramp;
        P_min  = physics.models{6}.P_min;
        if t_now < t_ramp
            P_now_for_phi = P_min + (P_max - P_min) * t_now / t_ramp;
        else
            % Use peak (not mean) cyclic pressure: peak drives crack opening
            P_now_for_phi = P_max;
        end
        % v11: no body_stress needed (degradation-driven model)
        
        % Solve
        solver.Solve();
        % ========== 实时监控 ==========
   if isfinite(verboseEvery) && verboseEvery > 0 && mod(tstep, verboseEvery) == 0
        [phiType, phiStep] = physics.dofSpace.getDofType({"phi"});
        phiNodes = mesh.GetAllNodesForGroup(1);
        phiDofs = physics.dofSpace.getDofIndices(phiType, phiNodes);
        max_phi = max(physics.StateVec{phiStep}(phiDofs));
        max_H   = max(physics.models{1}.Hist(:));
        CP_ratio = physics.models{1}.ContactPressureRatio;
        fprintf('[%s] Step %4d | Time %.1f days | φ_max=%.4f | H_max=%.2e | CP_ratio=%.3f\n', ...
            datestr(now,'HH:MM:SS'), tstep, physics.time/86400, max_phi, max_H, CP_ratio);
    end
    % ===================================================
        % Update time
        physics.time = physics.time + physics.dt;
        tvec(end+1) = tvec(end) + physics.dt;
        
        % Track results
        results.CL_avg(end+1)    = physics.models{2}.CL_int / mesh.Area(1);
        results.CL_max(end+1)    = physics.models{2}.CL_max;
        results.LFrac(end+1)     = physics.models{1}.LFrac;
        results.degradation_connected(end+1) = physics.models{1}.SpanningPath;
        results.alpha_avg(end+1) = physics.models{7}.alpha_avg;
        results.alpha_max(end+1) = physics.models{7}.alpha_max;
        results.totalCycles(end+1) = physics.models{8}.totalCycles;
        
        % Max fatigue damage from integration-point history avoids nodal
        % projection overshoot in the L2 recovery of D_f.
        if isprop(physics.models{8}, 'Df_max') && ~isempty(physics.models{8}.Df_max)
            results.Df_max(end+1) = min(1, physics.models{8}.Df_max);
        else
            allNodes = mesh.GetAllNodesForGroup(1);
            Df_dofs = physics.dofSpace.getDofIndices(...
                physics.dofSpace.getDofType({"D_f"}), allNodes);
            results.Df_max(end+1) = min(1, max(physics.StateVec{4}(Df_dofs)));
        end

        screening = collect_screening_metrics(physics, mesh);
        results.phi_avg(end+1) = screening.phi_avg;
        results.phi_p99(end+1) = screening.phi_p99;
        results.alpha_p99(end+1) = screening.alpha_p99;
        results.Df_p99(end+1) = screening.Df_p99;
        results.screening_index(end+1) = screening.index;
        
        if leakageModelIndex > 0
            leakage = physics.models{leakageModelIndex}.Evaluate(physics);
            results.leakage_Q(end+1) = leakage.Q;
            results.leakage_Q_poiseuille(end+1) = leakage.Q_poiseuille;
            results.leakage_kmax(end+1) = leakage.maxPermeability;
            results.leakage_aperture_max(end+1) = leakage.maxAperture;
            results.leakage_path_cost(end+1) = leakage.pathCost;
            results.leakage_path_mean_aperture(end+1) = leakage.pathMeanAperture;
            results.leakage_path_min_aperture(end+1) = leakage.pathMinAperture;
            results.leakage_connected(end+1) = leakage.connectedPath;
            results.contactPressureMean(end+1) = leakage.contactPressureMean;
            results.contactClosureMean(end+1) = leakage.contactClosureMean;
        end
        Q_allow = get_option(results, 'Q_allow', get_option(overrides, 'Q_allow', 1e-12));
        if isfield(results, 'leakage_Q') && ~results.leakage_failed_darcy(end) && ...
                results.leakage_Q(end) >= Q_allow
            results.leakage_failed_darcy(end+1) = true;
            results.leakage_time_darcy = physics.time / 86400;
        else
            results.leakage_failed_darcy(end+1) = results.leakage_failed_darcy(end);
        end
        if isfield(results, 'leakage_Q_poiseuille') && ~results.leakage_failed_poiseuille(end) && ...
                results.leakage_Q_poiseuille(end) >= Q_allow
            results.leakage_failed_poiseuille(end+1) = true;
            results.leakage_time_poiseuille = physics.time / 86400;
        else
            results.leakage_failed_poiseuille(end+1) = results.leakage_failed_poiseuille(end);
        end
        contactMetrics = collect_contact_metrics(physics);
        results.contact_top_mean(end+1) = contactMetrics.topMean;
        results.contact_bottom_mean(end+1) = contactMetrics.bottomMean;
        results.contact_right_mean(end+1) = contactMetrics.rightMean;
        results.contact_top_max(end+1) = contactMetrics.topMax;
        results.contact_bottom_max(end+1) = contactMetrics.bottomMax;
        results.contact_right_max(end+1) = contactMetrics.rightMax;
        
        % Print progress
        fprintf('  LFrac=%.3fmm  CL_avg=%.4f  alpha_max=%.4f  Df_max=%.4f', ...
            results.LFrac(end)*1e3, results.CL_avg(end), ...
            results.alpha_max(end), results.Df_max(end));
        if isfield(results, 'leakage_Q')
            fprintf('  Q_D=%.3e  Q_P=%.3e m^3/s  leakPath=%d  pcTop=%.2fMPa', ...
                results.leakage_Q(end), results.leakage_Q_poiseuille(end), ...
                logical(results.leakage_connected(end)), ...
                results.contact_top_mean(end)/1e6);
        end
        fprintf('\n');
        
        % Plot every 20 steps
        if isfinite(plotEvery) && plotEvery > 0 && mod(tstep, plotEvery) == 0
            plotres_seal(physics, tvec, results, mesh);
        end
        
        % Save every 50 steps
        if isfinite(saveEvery) && saveEvery > 0 && mod(tstep, saveEvery) == 0
            filename = savefolder + string(tstep);
            save(filename, "mesh","physics","solver","dt","tvec",...
                "results","n_max","tmax");
        end
        
        % Termination checks
        if physics.time > tmax
            disp("Reached maximum simulation time.");
            break;
        end
        
        % Do not terminate on the auxiliary connected-degradation extent.
        sealWidth = max(mesh.Nodes(:,1)) - min(mesh.Nodes(:,1));
        if results.LFrac(end) > sealWidth * 0.8
            fprintf('WARNING: left-connected degradation exceeds 80%% of the width at t=%.1f days.\n', ...
                physics.time/86400);
        end
    end
    
    % Final save
    filename = savefolder + "end";
    save(filename, "mesh","physics","solver","dt","tvec",...
        "results","n_max","tmax");
    
    toc(tmr)
end

function value = get_option(s, name, defaultValue)
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end

function metrics = collect_screening_metrics(physics, mesh)
    allNodes = mesh.GetAllNodesForGroup(1);
    [types, steps] = physics.dofSpace.getDofType({"phi", "alpha_a"});
    phiDofs = physics.dofSpace.getDofIndices(types(1), allNodes);
    alphaDofs = physics.dofSpace.getDofIndices(types(2), allNodes);
    phi = min(max(physics.StateVec{steps(1)}(phiDofs), 0), 1);
    alpha = min(max(physics.StateVec{steps(2)}(alphaDofs), 0), 1);
    Df = physics.models{8}.Df_ip(:);

    metrics.phi_avg = mean(phi);
    metrics.phi_p99 = prctile(phi, 99);
    metrics.alpha_p99 = prctile(alpha, 99);
    metrics.Df_p99 = prctile(Df, 99);
    metrics.index = max([metrics.phi_p99, metrics.alpha_p99, metrics.Df_p99]);
end

function out = merge_struct(out, in)
    names = fieldnames(in);
    for k = 1:numel(names)
        out.(names{k}) = in.(names{k});
    end
end

function idx = find_model(physics, modelName)
    idx = 0;
    for m = 1:length(physics.models)
        if isprop(physics.models{m}, 'myName') && physics.models{m}.myName == modelName
            idx = m;
            return;
        end
    end
end

function metrics = collect_contact_metrics(physics)
    metrics.topMean = 0;
    metrics.bottomMean = 0;
    metrics.rightMean = 0;
    metrics.topMax = 0;
    metrics.bottomMax = 0;
    metrics.rightMax = 0;
    
    for m = 1:length(physics.models)
        mdl = physics.models{m};
        if ~(isprop(mdl, 'myName') && mdl.myName == "ContactPenalty")
            continue;
        end
        groupName = "";
        if isprop(mdl, 'myGroup')
            groupName = mdl.myGroup;
        end
        switch groupName
            case "Top"
                metrics.topMean = mdl.meanPressure;
                metrics.topMax = mdl.maxPressure;
            case "Bottom"
                metrics.bottomMean = mdl.meanPressure;
                metrics.bottomMax = mdl.maxPressure;
            case "Right"
                metrics.rightMean = mdl.meanPressure;
                metrics.rightMax = mdl.maxPressure;
        end
    end
end


%% =============================================
%%  PLOTTING FUNCTION
%% =============================================
function plotres_seal(physics, tvec, results, mesh)
    figure(42)
    clf
    set(gcf, 'Position', [50 50 1400 900])
    
    % Domain extents (in meters)
    Lx_dom = max(mesh.Nodes(:,1));
    Ly_dom = max(mesh.Nodes(:,2));
    
    %% Phase-field (damage/crack pattern)
    subplot(3,4,1)
        physics.PlotNodal("phi", 0, "Internal");
        view(2);
        title('\phi (damage)')
        colorbar
        caxis([0 1])
        axis equal
        xlim([0 Lx_dom]); ylim([0 Ly_dom]);
        xlabel('x [m]'); ylabel('y [m]');
    
    %% Hydrostatic stress (2D view!)
    subplot(3,4,2)
        physics.PlotIP("sh", "Internal");
        view(2);
        title('\sigma_H [Pa]')
        colorbar
        axis equal
        xlim([0 Lx_dom]); ylim([0 Ly_dom]);
        xlabel('x [m]'); ylabel('y [m]');
    
    %% Fluid concentration
    subplot(3,4,3)
        physics.PlotNodal("CL", 0, "Internal");
        view(2);
        title('C_L (fluid)')
        colorbar
        caxis([0 1])
        axis equal
        xlim([0 Lx_dom]); ylim([0 Ly_dom]);
        xlabel('x [m]'); ylabel('y [m]');
    
    %% Aging variable
    subplot(3,4,4)
        physics.PlotNodal("alpha_a", 0, "Internal");
        view(2);
        title('\alpha_{age}')
        colorbar
        caxis([0 1])
        axis equal
        xlim([0 Lx_dom]); ylim([0 Ly_dom]);
        xlabel('x [m]'); ylabel('y [m]');
    
    %% Fatigue damage
    subplot(3,4,5)
        physics.PlotNodal("D_f", 0, "Internal");
        view(2);
        title('D_f (fatigue)')
        colorbar
        caxis([0 1])
        axis equal
        xlim([0 Lx_dom]); ylim([0 Ly_dom]);
        xlabel('x [m]'); ylabel('y [m]');
    
    %% Deformed shape with auto-scaled displacement
    subplot(3,4,6)
        [~, dxSteps_u] = physics.dofSpace.getDofType({"dx","dy"});
        max_disp = max([max(abs(physics.StateVec{dxSteps_u(1)})), ...
                        max(abs(physics.StateVec{dxSteps_u(2)}))]);
        if max_disp > 1e-15
            dscale = round(0.05 * Lx_dom / max_disp);
            dscale = max(1, min(dscale, 5000));
        else
            dscale = 1;
        end
        physics.PlotNodal("phi", dscale, "Internal");
        view(2);
        title(['\phi (deformed x' num2str(dscale) ')'])
        colorbar
        caxis([0 1])
        axis equal
        xlabel('x [m]'); ylabel('y [m]');
    
    %% Time histories
    t_days = tvec / (3600*24);
    
    subplot(3,4,7)
        plot(t_days, results.CL_avg, 'b-', 'LineWidth', 1.5)
        xlabel('Time [days]')
        ylabel('Avg C_L [-]')
        title('Fluid Uptake')
        grid on
    
    subplot(3,4,8)
        plot(t_days, results.LFrac*1000, 'r-', 'LineWidth', 1.5)
        xlabel('Time [days]')
        ylabel('L_{frac} [mm]')
        title('Crack Length')
        grid on
    
    subplot(3,4,9)
        plot(t_days, results.alpha_avg, 'g-', 'LineWidth', 1.5)
        hold on
        plot(t_days, results.alpha_max, 'g--', 'LineWidth', 1)
        xlabel('Time [days]')
        ylabel('\alpha_{age} [-]')
        title('Aging')
        legend('Average','Maximum','Location','best')
        grid on
    
    subplot(3,4,10)
        plot(t_days, results.Df_max, 'm-', 'LineWidth', 1.5)
        xlabel('Time [days]')
        ylabel('D_f max [-]')
        title('Fatigue Damage')
        grid on
    
    subplot(3,4,11)
        if isfield(results, 'leakage_Q')
            semilogy(t_days, max(results.leakage_Q, realmin), 'k-', 'LineWidth', 1.5)
            hold on
            leg = ["Darcy"];
            if isfield(results, 'leakage_Q_poiseuille')
                semilogy(t_days, max(results.leakage_Q_poiseuille, realmin), 'k--', 'LineWidth', 1.1)
                leg(end+1) = "Poiseuille";
            end
            if isfield(results, 'Q_allow')
                yline(results.Q_allow, 'r:', 'Q_{allow}', 'LineWidth', 1.0)
                leg(end+1) = "Q_{allow}";
            end
            xlabel('Time [days]')
            ylabel('Q [m^3/s]')
            title('Leakage Rate')
            legend(leg, 'Location', 'best')
        else
            plot(results.totalCycles, results.LFrac*1000, 'k-', 'LineWidth', 1.5)
            xlabel('Cycles [-]')
            ylabel('L_{frac} [mm]')
            title('Crack vs Cycles')
        end
        grid on
    
    %% Screening summary
    subplot(3,4,12)
        if isfield(results, 'contact_top_mean')
            plot(t_days, results.contact_top_mean/1e6, 'Color', [0.85 0.25 0.15], 'LineWidth', 1.4)
            hold on
            plot(t_days, results.contact_bottom_mean/1e6, 'Color', [0.10 0.45 0.75], 'LineWidth', 1.2)
            plot(t_days, results.contact_right_mean/1e6, 'Color', [0.20 0.60 0.30], 'LineWidth', 1.2)
            xlabel('Time [days]')
            ylabel('Mean p_c [MPa]')
            title('Contact Pressure')
            legend('Top','Bottom','Right','Location','best')
            grid on
        else
            screening_index = results.screening_index(end);
            bar_data = [results.phi_p99(end)*100, ...
                        results.Df_p99(end)*100, ...
                        results.alpha_p99(end)*100, ...
                        results.CL_avg(end)*100];
            bar(bar_data)
            set(gca, 'XTickLabel', {'phi P99', 'Fatigue P99', 'Aging P99', 'Fluid avg'})
            ylabel('Degradation [%]')
            title(['Screening index (I_S=' sprintf('%.2f', screening_index) ')'])
            ylim([0 100])
            grid on
        end
    
    sgtitle(['Seal Analysis - t = ' ...
        sprintf('%.2f', tvec(end)/3600/24) ' days, ' ...
        sprintf('%.0f', results.totalCycles(end)) ' cycles'])
    
    drawnow();
end
