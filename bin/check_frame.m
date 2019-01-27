% ****************==- �����ļ��� -==*********************
if length(result_folder)==0
   error('  *** please give one filefolder name for result ***'); 
end

if ~(result_folder(end) == '\' || result_folder(end) == '/')
   error('  *** please end the result_folder with / or \ ***'); 
end
if ~exist(result_folder,'dir')
    mkdir(result_folder);
    copyfile('doc/Result.xls',result_folder)
    copyfile('doc/Tex.xls',result_folder)
end

% ****************==- ������ -==*********************
if length(func_name)==0
   error('  *** please give one test function name  ***'); 
end

% ****************==- �㷨�� -==*********************
if length(algo_name)==0
   error('  *** please give one Algorithm name ***'); 
end

if  ~exist( [result_folder 'boxfig/'] ,'dir') mkdir([result_folder 'boxfig/']);    end
if  ~exist( [result_folder 'data/'] ,'dir')   mkdir([result_folder 'data/']);      end
