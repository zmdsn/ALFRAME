function [ algo_config , func_config, objfunc, PopuSize,maxFES,DIM,Xmin,Xmax , f, init] = get_parameters( func_name , algo_config,func_config,default_set )
% This function is used to preprocess the parameters required by 
% the algorithm

% output :
%   algo_config    struct
%   PopuSize
%   maxFES 

    global func_name_flag
    if ~isequal(func_name_flag,func_name)
        % set PopuSize
        if isfield(algo_config,'PopuSize')
            PopuSize = algo_config.PopuSize;
        else
            PopuSize = default_set.PopuSize;
            algo_config.PopuSize = PopuSize;
        end

        % set maxFES
        if isfield(algo_config,'maxFES')
            maxFES = algo_config.maxFES;
        else
            maxFES = default_set.maxFES;
            algo_config.maxFES = maxFES;
        end
        
        % set objfunc
        if isa(func_name,'function_handle') 
            func_config.handle = func_name;
            objfunc = feval('func_handle',func_config);  
        else
            objfunc = feval(func_name,func_config);  
        end
        
        DIM  = objfunc.Xdim;          
        Xmin = repmat(objfunc.Xmin,1,PopuSize);
        Xmax = repmat(objfunc.Xmax,1,PopuSize);
        f = @objfunc.obj_value;
        init = @objfunc.initialize;

    end
end

