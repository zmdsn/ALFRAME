AlgoFrame使用手册
---


[TOC]


### 前言和简介
AlgoFrame,优化算法开发框架.在这里将写明为什么要开发AlgoFrame,以及开发怎样的AlgoFrame,希望读者能够认真的阅读,这样会让你省去大量的试错时间.
开发AlgoFrame最初的目的是可以方便的进行群智能算法的测试和开发.不过渐渐的发现其实可以用于不同的情境中,比如传统优化算法的测试和开发、分类聚类、图像处理等等的情境中.能够胜任这些任务最主要的原因是开发框架使用了matlab中的类。由于在matlab2015b版本中，对类的使用进行了优化，使之运行速度加快，所以最好要使用matlab2015b及以上的版本。对matlab类不需要太高深的理解就能使用，这是一大优点。如果对matlab类不了解可以参见附录中关于类的简单介绍。
一个新算法的开发往往要经过以下几个内容：
> 
* 构思算法蓝图--即对算法进行一个全局的构思，提出算法的大致轮廓
* 实现算法观察算法的结果，并且调整构思蓝图
* 查看算法的收敛曲线，并且调整构思蓝图
* 将算法与别的算法进行比较
* 算法参数的微调,和讨论
* 使用算法求解不同维度的问题以查看算法对低维和高维不同问题的性能

针对以上的不同情境,我们希望都能有一个快速的解决方案.
#### 算法的不同阶段的解决方案
##### 算法开发和调试阶段
在实现一个算法的时候需要调试算法,在调整某些部分之后常常希望能够快速的知道调整之后的算法是否有可行这个时候你可以使用如下的命令:
```matlab
PSO('f1')
```
得到的结果就是粒子群算法对问题f1的最优函数值。在调试算法的时候使用这个命令可以很快地知道自己错在了哪里,并及时的进行修改.由于我们寻找的是函数的最小值，所以这样可以很方便的判断出算法是否达到了理想中的效果。这种使用算法名字和函数名字调用的方式是该框架最大的特点。使你在测试阶段能够快速的完成算法的开发.
##### 查看算法的收敛曲线
在算法开发完成之后需要查看这个算法的收敛曲线,并通过曲线观察算法在不同阶段的表现,并针对做出调整,这时我们可以使用以下的命令:
```matlab
altest('PSO','f1',1);
```
![](img\曲线.png)


##### 将算法与别的算法进行比较
在不同算法比较的时候就需要修改脚本文件comp.m中的如下部分:
```matlab
% ****************==- 框架参数设置 -==*********************
Run = 30;                   % 每个算法单独运行次数
result_folder = 'test2\';   % 结果文件保存目录 以\或/ 结尾

% ****************==- 算法参数设置 -==*********************
algo_name = { 'PSO','bRGA','DE1'};  % 算法名称
algo_config.PopuSize = 50;          % 种群规模
algo_config.maxFES = 100000;        % 目标函数运行次数
nAlgorithm = length(algo_name);

% ****************==- 测试函数设置 -==*********************
% func_name = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'};
									% 需要求解的问题名称
func_config.Xdim = 30;				% 需要求解的问题的维度
nFunc_num = length(func_name);
```

在设置的文件夹下会得到多个测试函数的对比图像和两个表格文件.
Result.xls用于记录运行最终数据如下.

| min | mean | median | std | t-test |
|--------|--------|--------|--------|
|1|1|1|1|1|
|0|0|0|0|0|

Tex.xls用于记录测试函数的表达式

| 表达式 |
|--------|
|y = x^2|
|y = x^3|

两个文件的表头均未给出,需要自行添加


##### 其它
在使用算法参数的微调,查看算法求解不同维度的性能的时候需要算法本身做一些调整,这本分内容将在讲完如何开发新算法和添加函数之后再讲述.

#### 框架的可扩展性
算法开发框架目前看似庞大,不过基本的内容可以很方便的拆分和组装,可以很方便的添加和移除某些内容.算法的目录和基本的功能如下
> —— algo
> ———— bin  		% 算法框架相关函数
> ———— doc  		% 算法相关的文档
> ———— algorithms 	% 存放各种算法的目录
> ———— functions  	% 存放各种测试函数的目录
> —————— @func_class  % 测试函数的基类
> —————— cec14 等文件夹   % cec测试函数 
> ———— tools      	% 一些有用的小工具
> ———— GIFS      	% 存放gif动态图


### 添加一个新的测试函数
添加一个新的测试函数可以将doc文件夹中的func_name.m文件中的相关设置修改就可以了.其中需要将func_name修改为你想要的名字并且使用该名字文件名.
```matlab
classdef func_name < func_class 
   methods
        %% 主运算函数
        function [f, g, h] = object_function(obj, x , specific )
            % 注意这里x为列向量
            % specific 预留值,建议使用结构数组,用于传递额外的参数,一般不用
            % f(x) = 0
            % g(x) <= 0
            % h(x) = 0
            g = zeros(obj.Fneq,1);  % 这个地方是不等式约束
            h = zeros(obj.Feq,1);	% 这个地方是等式约束	
			f = sum(x.^2) ;			% 这个地方是函数的具体表达式
        end
        
       %% 构造函数 constructed function
        function obj = func_name(specific)
            % specific 预留值,使用结构数组,用于传递额外的参数
            % 设置函数参数的默认值
            if ~exist('specific','var')
                specific = [];  % 参数初始化
            end
            if ~isfield(specific,'Xdim') specific.Xdim = 30;  end % 变量维数
            if ~isfield(specific,'Xmin') specific.Xmin = -100*ones(specific.Xdim,1); end  % 下界 
            if ~isfield(specific,'Xmax') specific.Xmax = 100*ones(specific.Xdim,1);  end   % 上界
            if ~isfield(specific,'Fmin') specific.Fmin = 0;  end  % 目标函数最小值
            if ~isfield(specific,'Tex')  specific.Tex  = 'f=\sum\limits_{i=1}^D x_i^2';  end   % 函数tex表达式
            if ~isfield(specific,'Nobj') specific.Nobj = 1;  end  % 目标函数值个数
            if ~isfield(specific,'Fneq') specific.Fneq = 0;  end  % 不等式约束g
            if ~isfield(specific,'Feq')  specific.Feq  = 0;  end  % 等式约束h
            obj = obj@func_class(specific);
        end
   end
end
```


将以上各个部分修改完成之后放到相应的目录下即可.
#### 测试函数相关的应用
##### 调用一个测试函数
由于测试函数时一个类,调用一个测试函数也就是使用初始化一个类,使用如下的方式即可完成类的初始化.
```matlab
objfunc = feval('f1');
或者
objfunc = f1;
```
##### 绘制测试函数的性质曲线
fun_figure(x)
1 维曲线
2 二维曲面和等高线
3 等高线
例如运行以下命令就会得到下边的图像:
```matlab
objfunc.fun_figure(2) ;
```

![](img\func.png)

##### 运行时的种群动态图
将函数插入到算法中,算法就会以2维显示在运行时种群的动态图像.
```matlab
objfunc.draw_running(swarm,'PSO');	  % 会保gif格式图片到GIFS文件夹下
objfunc.draw_running(swarm);			% 不会保存
```
![](img\PSO.gif)

##### 初始化种群
将使用公式$X^i = rand*(UB - LB) + LB$ 的方式产生初始种群
```matlab
objfunc.initialize( PopuSize )
```
这个命令会产生PopuSize个个体组成的一个矩阵,每一列为一个解,一共PopuSize列,即PopuSize个解.



### 添加一个新的算法
有的时候算法是需要自己写,有的时候可能是从某个地方获取到的,无论怎样,把算法添加到框架是一个必要的步骤.首先将算法框架传入参数转化为适应于具体算法的参数.如下是PSO的一种转换方式,首先将函数的定义调整如下,

```matlab
function [best_fvalue,best_solution,run_series,run_info,user_set] = PSO(func_name , algo_config,func_config,user_config)
```
其中输入:
* func_name   (字符串) 是需要测试函数名
* algo_config (结构体) 是算法设置
* func_config (结构体) 是测试函数设置
* user_config (结构体,一般不用) 是用户自定义数据.
输出:
* best_fvalue   (数值) 全局最优函数值
* best_solution (向量) 全局最优解
* run_series    (m*2的矩阵,第一列记录测试函数计算次数,第二列记录目前算法得到的最优解) 算法迭代过程中记录的每代最优解(用于绘制算法收敛曲线)
* run_info  	运行信息
* user_set		用户自定义输出(一般不用)


其次,对于一般的代码,需要在代码前使用如下的初始化设置取代原来算法的初始化设置,将迭代次数和种群规模等转换为已有算法的设置.如果是新增一个算法,复制过去,改一下算法名称就可以了.
```matlab
function [best_fvalue,best_solution,run_series,run_info,user_set] = PSO(func_name , algo_config,func_config,user_config)
    mytime = cputime;
    format long;
    format compact;
    %% *****************==- 初始化测试函数 -==**********************
    if exist('func_config','var')
        objfunc = feval(func_name,func_config);  % 实例化一个目标函数类
    else
        objfunc = feval(func_name);  % 实例化一个目标函数类
        func_config = [];
    end
    objfunc = feval(func_name,func_config);
    DIM = objfunc.Xdim;          % 获取目标函数类的一些参数
    Xmin =  objfunc.Xmin(1);
    Xmax = objfunc.Xmax(1);

    %% *****************==- 初始化算法参数 -==**********************
    if ~exist('algo_config','var')
        algo_config = [];  % 参数初始化
    end
    % 种群规模
    if isfield(algo_config,'PopuSize')
        popsize = algo_config.PopuSize;
    else
        popsize = 10;
    end
    % 迭代次数
    if isfield(algo_config,'maxFES')
        maxFES = algo_config.maxFES;
    else
        maxFES = 60000;
    end

```

再次,在算法的最后进行算法结果数据的调整,使得之符合数据输出的格式
```matlab
    %% ****************==- 算法结果数据处理 -==*********************
    [best_fvalue,minidx]= min(fit);
    best_solution = uSet(minidx, :)';
    run_series = outcome;
    Altime = cputime - mytime;                  % 计算时间
    run_info = [Altime];

```

最后,将算法内的用于求解函数值的地方统一使用如下代码代替,其中输入要求为**列向量**或者是一个矩阵,**每一列代表一个解**.结果输出为一个行向量,对应的为每个解的目标函数值.
```matlab
y = objfunc.obj_value(X);
```
注意原算法得到结果是行向量还是列向量,如果行和列向量不一致会导致错误.


### 其它的一些说明
#### 关于CEC系列测试函数的说明
CEC系列的测试函数是源于CEC会议,具体的代码可以到学者Suganthan的网站下载http://www.ntu.edu.sg/home/epnsugan/
这一套测试函数是这个领域里边比较权威的一套测试函数,包含单目标,多目标,约束和动态等各种类型的问题,每年会举行一次算法的竞赛,用于推动智能算法的研究.
现在对需要使用的技术进行一些简单的说明.
需要使用到matlab和c语言的混合编程,不过不用担心,挺简单的,主要内容就是学会matlab中mex的运用.在下载的文件夹中会有readme.txt文件,里边有具体的使用方法.如果需要详细学习可以到http://blog.sina.com.cn/s/blog_468651400100coas.html去具体学习一下.不过在Suganthan的网站上可以下载好matlab版本的,具体的机器上需要重新编译一下,所以此处才会提到这个问题.
以CEC14为例,在代码中我修改了文件的路径,更改了cec14的input文件路径,使之适应我们的算法框架,不过这没有什么关系.在matlab2015b之前的某一版本开始,就不需要安装mex了,这大大简化了我们的操作.
在matlab命令行中输入一下命令 
首先使用命令 mex -setup 选择编译器
之后转到cec14_func.cpp所在的文件夹使用  mex cec14_func.cpp 编译完成得到一个类似cec14_func.mexw64的文件即可

有的时候可能提示你没有安装任何编译器,这时,你需要安装一个C语言的编译器.然后再选择编译器

**说明**:
cec系列问题可以直接调用,并不一定需要使用函数类的继承,具体请参照各个不同测试函数的readme.txt

#### 关于TSP问题的测试函数
算法框架里放置了常用的TSP问题测试数据,并简单的写了一个使用这些数据的一个测试函数,如果需要的话可以自行拓展.
在测试函数里也添加了相应的函数支持,例如:
```matlab
    % draw_road(obj,citys,road,specific)            绘制路径图(TSP)
    % is_one2n(obj,x)                               判断是否为1:n的数列的乱序
```


#### 关于多目标和带约束问题的算法
测试函数已经给出,具体的算法还没有完成,暂时还没有完工.




### 离线测试框架
在使用离线测试的时候首先每个算法每个函数独立运行60次，每次在60次结果中选择所需要次数的结果作为该次运行结果。

f 运行中必须添加func_config





