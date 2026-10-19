function opts = set_opts_fmincon(n_elements)

opts = optimoptions('fmincon');

% shared options
opts.Algorithm = 'interior-point';
%opts.CheckGradients=true;
opts.ConstraintTolerance = 1e-4;
opts.Diagnostics = 'on';
opts.Display='iter-detailed';
opts.FiniteDifferenceType = 'central';
opts.MaxIterations = 600;
opts.MaxFunctionEvaluations = 3001;
opts.OptimalityTolerance = 1e-5; 
opts.SpecifyConstraintGradient = true;
opts.SpecifyObjectiveGradient = true;
opts.StepTolerance = 1e-6;
opts.TypicalX =[5*ones(n_elements,1);1];

% IP options    
opts.BarrierParamUpdate = 'predictor-corrector'; % 'monotone'
opts.HessianApproximation = 'lbfgs'; 
opts.HonorBounds = true;
opts.InitBarrierParam = 1e-3; 
opts.InitTrustRegionRadius = 1e4;
opts.SubproblemAlgorithm = 'cg';
% opts.EnableFeasibilityMode = true;