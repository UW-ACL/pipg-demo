function [] = codegen_epipg(pipg_dvec_pbm,flag,varargin)

    if nargin <= 2
        mxcfg = coder.config('mex');
        mxcfg.DynamicMemoryAllocation                       = 'Off';
        mxcfg.EnableVariableSizing                          = false;
        mxcfg.PreserveArrayDimensions                       = false;
        mxcfg.ResponsivenessChecks                          = false;
        mxcfg.RowMajor                                      = true;
        mxcfg.SaturateOnIntegerOverflow                     = false;
        mxcfg.IntegrityChecks                               = false;
        mxcfg.HighlightPotentialRowMajorIssues              = false;
        mxcfg.HighlightImplicitExpansionIssues              = false;
        mxcfg.GlobalDataSyncMethod                          = 'NoSync';
        mxcfg.ConstantInputs                                = 'IgnoreValues';
        mxcfg.TargetLang                                    = 'C';

        mxcfg.EnableAutoParallelization                     = false;
        mxcfg.EnableAutoExtrinsicCalls                      = false;
        mxcfg.ExtrinsicCalls                                = false; 
        mxcfg.EnableRuntimeRecursion                        = false; %%%
        mxcfg.CompileTimeRecursionLimit                     = 0; %%
    else
        mxcfg = varargin{1};
    end
    
    pipg_pbm = pipg_dvec_pbm;
    pipg_args = {pipg_pbm.qinit, pipg_pbm.uinit, pipg_pbm.vinit, pipg_pbm.err,...
                 coder.Constant(pipg_pbm.eps_pow), coder.Constant(pipg_pbm.eps_pipg), coder.Constant(pipg_pbm.eps_infeas), coder.Constant((pipg_pbm.max_iter_pow)), coder.Constant((pipg_pbm.max_iter)), coder.Constant(pipg_pbm.test_freq), coder.Constant(pipg_pbm.rho), (pipg_pbm.omg), coder.Constant(pipg_pbm.lam),...
                 ((pipg_pbm.N)), ((pipg_pbm.M)), coder.Constant(pipg_pbm.Q), coder.Constant(pipg_pbm.R), coder.Constant(pipg_pbm.qq), coder.Constant(pipg_pbm.ru), coder.Constant(pipg_pbm.A), coder.Constant(pipg_pbm.B), coder.Constant(pipg_pbm.AT), coder.Constant(pipg_pbm.BT),...
                 coder.Constant(pipg_pbm.q0), coder.Constant(pipg_pbm.qf), coder.Constant(pipg_pbm.qmax), coder.Constant(pipg_pbm.umax)...
                 , true...
                 };

    switch flag
        case 0 % only epipg
            pipg_args(4) = [];
            pipg_args(end) = [];
            codegen epipg -config mxcfg -d solvers/epipg/codegen/mex/epipg -o solvers/epipg/epipg_mex -args pipg_args    
        case 1 % only epipg_verb
            codegen epipg_verb -config mxcfg -d solvers/epipg/codegen/mex/epipg_verb -o solvers/epipg/epipg_verb_mex -args pipg_args   
        case 2 % all
            codegen epipg_verb -config mxcfg -d solvers/epipg/codegen/mex/epipg_verb -o solvers/epipg/epipg_verb_mex -args pipg_args
    
            pipg_args(4) = [];
            pipg_args(end) = [];
            codegen epipg -config mxcfg -d solvers/epipg/codegen/mex/epipg -o solvers/epipg/epipg_mex -args pipg_args            
        otherwise
            error('Invalid flag');
    end

end