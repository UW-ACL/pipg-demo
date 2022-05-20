function [] = codegen_expipg_constsize(pp,ppv,alg_type)
%{
05/05/2022
Purnanand Elango

Generate C code for exPIPG functions written scripts MATLAB using the Coder
Variable-size data is *not* allowed; code generation is required everytime
problem size changes

Creates functions with suffixes of the form 
..._n?_m?_N?_mex

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
%}

    mxcfg = coder.config('mex');
    mxcfg.EnableVariableSizing                          = false;
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

    n = pp.n;                                           % State dim
    m = pp.m;                                           % Control dim
    N = pp.N;                                           % Horizon length
    Z = (n+m)*N+n;                                      % Length of vectorized decision vector
    Y = n*N;

    % Segment of the file name which will denote fixed problem size for
    % which the generated file is valid
    name_suffix = "_n" + num2str(n) + "_m" + num2str(m) + "_N" + num2str(N) + "_mex";

    switch alg_type{1}
        case 'VEC'
            % Creation of sparse objects requires (variable-size data) dynamic memory allocation
            mxcfg.EnableVariableSizing                          = true;
            mxcfg.DynamicMemoryAllocation                       = 'Threshold';
            mxcfg.DynamicMemoryAllocationForVariableSizeArrays  = 'Threshold';
            % Estimate static allocation threshold from memory occupied by matrix H
            H = ppv.H;
            dataH = whos('H');
            mxcfg.DynamicMemoryAllocationThreshold = 5*dataH.bytes;

            vec_args =  {coder.typeof(0.0,[Z 1],[0 0]),coder.typeof(0.0,[Y 1],[0 0]),...
                         coder.typeof(0.0,[Z 1],[0 0]),coder.typeof(0.0,[Y 1],[0 0]),...
                         coder.typeof(0.0,[Z 1],[0 0]),coder.typeof(0.0,[Y 1],[0 0]),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(0,1,0),coder.typeof(0,1,0),coder.typeof(0,1,0),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(ppv.P,[Z Z],[0 0]),coder.typeof(ppv.H,[Y Z],[0 0]),coder.typeof(ppv.HT,[Z Y],[0 0]),...
                         coder.typeof(0.0,[Z 1],[0 0]),coder.typeof(0.0,[Z 1],[0 0])};

            switch alg_type{2}
                case 'old_infeas'
                    func_name = "solvers/expipg/expipg_vec_v2" + name_suffix;
                    codegen("solvers/expipg/expipg_vec_v2.m","-config",mxcfg,"-d","solvers/expipg/codegen/mex/expipg_vec_v2/","-o",func_name,"-args",vec_args);
                case 'new_infeas'
                    func_name = "solvers/expipg/expipg_vec" + name_suffix;
                    codegen("solvers/expipg/expipg_vec.m","-config",mxcfg,"-d","solvers/expipg/codegen/mex/expipg_vec/","-o",func_name,"-args",vec_args);
                otherwise
                    error("Invalid infeasibility detection critiera; should be 'old_infeas' or 'new_infeas'");
            end
        case 'DVEC'

            dvec_args = {coder.typeof(0.0,[n N+1],[0 0]),coder.typeof(0.0,[m N],[0 0]),coder.typeof(0.0,[n N],[0 0]),...
                         coder.typeof(0.0,[n N+1],[0 0]),coder.typeof(0.0,[m N],[0 0]),coder.typeof(0.0,[n N],[0 0]),...
                         coder.typeof(0.0,[n N+1],[0 0]),coder.typeof(0.0,[m N],[0 0]),coder.typeof(0.0,[n N],[0 0]),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0,1,0),coder.typeof(0,1,0),coder.typeof(0,1,0),...
                         coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),...
                         coder.typeof(0,1,0),coder.typeof(0.0,[n n],[0 0]),coder.typeof(0.0,[m m],[0 0]),coder.typeof(0.0,[n n],[0 0]),coder.typeof(0.0,[n n],[0 0]),...
                         coder.typeof(0.0,[n m],[0 0]),coder.typeof(0.0,[m n],[0 0]),...
                         coder.typeof(0.0,[n 1],[0 0]),coder.typeof(0.0,[n 1],[0 0]),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0),coder.typeof(0.0,1,0)};
            
            switch alg_type{2}
                case 'old_infeas'
                    func_name = "solvers/expipg/expipg_dvec_v2" + name_suffix;
                    codegen("solvers/expipg/expipg_dvec_v2.m","-config",mxcfg,"-d","solvers/expipg/codegen/mex/expipg_dvec_v2/","-o",func_name,"-args",dvec_args);
                case 'new_infeas'
                    func_name = "solvers/expipg/expipg_dvec" + name_suffix;
                    codegen("solvers/expipg/expipg_dvec.m","-config",mxcfg,"-d","solvers/expipg/codegen/mex/expipg_dvec/","-o",func_name,"-args",dvec_args);
                otherwise
                    error("Invalid infeasibility detection critiera; should be 'old_infeas' or 'new_infeas'");
            end
        otherwise
            error("Incorrect algorithm type; should be 'VEC' or 'DVEC'");
    end

    % Remove codegen source file; only retain the binaries
    rmdir solvers/expipg/codegen/ s
end