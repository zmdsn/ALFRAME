for i = 2:15
    filename = ['D:\workspace\algo2.0\Test\OLPSO\cec15_f' num2str(i) '_30d_50.mat'];
    load(filename);
    run_data = run_data(1:end-1,:);
    save( filename , 'run_data');
end