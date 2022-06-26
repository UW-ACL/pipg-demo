function [] = codegen_expipg_varsize(pp,ppv,alg_type)
%{
05/05/2022
Purnanand Elango

Generate C code for exPIPG functions written scripts MATLAB using the Coder
Variable-size data is allowed so that a single generated code can handle
multiple problem sizes

Refer to the following pages for information about handling variable size input:
    - https://www.mathworks.com/help/coder/ug/code-generation-for-variable-size-data.html
    - https://www.mathworks.com/help/coder/ug/what-is-variable-size-data.html
    - https://www.mathworks.com/help/coder/ug/controlling-memory-allocation-of-bounded-data.html
    - https://www.mathworks.com/help/coder/ug/variable-size-data-when-dynamic-memory-allocation-is-disabled.html
    - https://www.mathworks.com/help/coder/ref/coder.typeof.html

Input:
    alg_type{1}:
        Vectorized ('VEC')
        De-vectorized ('DVEC')
    alg_type{2}:
        Infeasibility detection based on dual error ('old_infeas')
        New criteria developed by Yue Yu on 04/30/2022 ('new_infeas')
        New infeasibility and feasibility criteria developed by Yue Yu on 04/30/2022 and 06/13/2022 respectively ('new')
%}

    mxcfg = coder.config('mex');
    mxcfg.EnableVariableSizing                          = true;
    mxcfg.DynamicMemoryAllocation                       = 'Off';
    mxcfg.DynamicMemoryAllocationForVariableSizeArrays  = 'Never';
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
    mxcfg.EnableRuntimeRecursion                        = false;
    mxcfg.CompileTimeRecursionLimit                     = 0;

    n_max = pp.n_max;                                   % Max state dim
    m_max = pp.m_max;                                   % Max control dim
    N_max = pp.N_max;                                   % Max horizon length
    Z_max = (n_max+m_max)*N_max+n_max;                  % Max length of vectorized decision vector
    Y_max = n_max*N_max;

    switch alg_type{1}
        case 'VEC'
            % Creation of sparse objects requires (variable-size data) dynamic memory allocation
            mxcfg.DynamicMemoryAllocation                       = 'Threshold';
            mxcfg.DynamicMemoryAllocationForVariableSizeArrays  = 'Threshold';
            % Estimate static allocation threshold from memory occupied by matrix H
            H = ppv.H;
            dataH = whos('H');
            mxcfg.DynamicMemoryAllocationThreshold = 5*dataH.bytes;

            vec_args =  {coder.typeof(0.0,[Z_max 1],[1 0]),coder.typeof(0.0,[Y_max 1],[1 0]),...
                         coder.typeof(0.0,[Z_max 1],[1 0]),coder.typeof(0.0,[Y_max 1],[1 0]),...
                         coder.typeof(0.0,[Z_max 1],[1 0]),coder.typeof(0.0,[Y_max 1],[1 0]),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(0,1,0),coder.typeof(0,1,0),coder.typeof(0,1,0),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(ppv.P,[Z_max Z_max],[1 1]),coder.typeof(ppv.H,[Y_max Z_max],[1 1]),coder.typeof(ppv.HT,[Z_max Y_max],[1 1]),...
                         coder.typeof(0.0,[Z_max 1],[1 0]),coder.typeof(0.0,[Z_max 1],[1 0])};

            switch alg_type{2}
                case 'old_infeas'
                    codegen solvers/expipg/expipg_vec_v2 -config mxcfg -d solvers/expipg/codegen/mex/expipg_vec_v2/ -o solvers/expipg/expipg_vec_v2_mex -args vec_args
                case 'new_infeas'
                    codegen solvers/expipg/expipg_vec -config mxcfg -d solvers/expipg/codegen/mex/expipg_vec/ -o solvers/expipg/expipg_vec_mex -args vec_args
                case 'new'
                    codegen solvers/expipg/expipg_vec_v3 -config mxcfg -d solvers/expipg/codegen/mex/expipg_vec_v3/ -o solvers/expipg/expipg_vec_v3_mex -args vec_args    
                otherwise
                    error("Invalid infeasibility detection critiera; should be 'old_infeas', 'new_infeas' or 'new'");
            end
        case 'DVEC'

            dvec_args = {coder.typeof(0.0,[n_max N_max+1],[1 1]),coder.typeof(0.0,[m_max N_max],[1 1]),coder.typeof(0.0,[n_max N_max],[1 1]),...
                         coder.typeof(0.0,[n_max N_max+1],[1 1]),coder.typeof(0.0,[m_max N_max],[1 1]),coder.typeof(0.0,[n_max N_max],[1 1]),...
                         coder.typeof(0.0,[n_max N_max+1],[1 1]),coder.typeof(0.0,[m_max N_max],[1 1]),coder.typeof(0.0,[n_max N_max],[1 1]),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0,1,0),coder.typeof(0,1,0),coder.typeof(0,1,0),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(0,1,0),coder.typeof(0.0,[n_max n_max],[1 1]),coder.typeof(0.0,[m_max m_max],[1 1]),coder.typeof(0.0,[n_max n_max],[1 1]),coder.typeof(0.0,[n_max n_max],[1 1]),...
                         coder.typeof(0.0,[n_max m_max],[1 1]),coder.typeof(0.0,[m_max n_max],[1 1]),...
                         coder.typeof(0.0,[n_max 1],[1 0]),coder.typeof(0.0,[n_max 1],[1 0]),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0)};
            
            switch alg_type{2}
                case 'old_infeas'
                    codegen solvers/expipg/expipg_dvec_v2 -config mxcfg -d solvers/expipg/codegen/mex/expipg_dvec_v2/ -o solvers/expipg/expipg_dvec_v2_mex -args dvec_args
                case 'new_infeas'
                    codegen solvers/expipg/expipg_dvec -config mxcfg -d solvers/expipg/codegen/mex/expipg_dvec/ -o solvers/expipg/expipg_dvec_mex -args dvec_args
                case 'new'
                    codegen solvers/expipg/expipg_dvec_v3 -config mxcfg -d solvers/expipg/codegen/mex/expipg_dvec_v3/ -o solvers/expipg/expipg_dvec_v3_mex -args dvec_args    
                otherwise
                    error("Invalid infeasibility detection critiera; should be 'old_infeas', 'new_infeas' or 'new'");
            end
        otherwise
            error("Incorrect algorithm type; should be 'VEC' or 'DVEC'");
    end
    
    % Remove codegen source file; only retain the binaries
    rmdir solvers/expipg/codegen/ s
end