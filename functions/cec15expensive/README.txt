1. Intro

This software package contain the Matlab wrapper for the C version code for 
"Problem Definitions and Evaluation Criteria for CEC 2015 Special Session 
on Bound Constrained Single-Objective Numerical Optimization Competition 
Part B: Computationally Expensive Single Objective Optimization".

For detail of the problems, please check out the related technical report.

2. Usage

The main interface of the software package defined in cec15problems.c. 
To use the functions, you have to compile your file with cec15_test_func.c 
and record_result_iml.c and rp_iml.c. 

For Windows platform, you have to run:

>> mex cec15problems.c cec15_test_func.c record_result_iml.c  -D_WINDOWS

For Linux-like platform, you have to run:

>> mex cec15problems.c cec15_test_func.c record_result_iml.c

If you don't want to use the result recording and statistics functionalities, 
just define NO_RECORDING when you compile cec15problems.c.

mex cec15problems.c cec15_test_func.c  -D_WINDOWS -DNO_RECORDING

3. Interface

The main part of them are three functions as follows:

/*
 * evaluate specific function
 */
cec15problems('eval',x,func_num)

x will contain the values need to be evaluated, the length of x will be taken
as the dimension of the problem, can be 10 and 30 here; 
func_num defines the function you want to evaluate, can be 1,2,3,..., 15.

/* Recording and statistics functions */
cec15problems('runnumber',runtime);
cec15problems('output',dirname, prefix);

These two functions can be turned of by define NO_RECORDING. 


All the evaluations after call cec15problems('runnumber',n) and before 
the next call of it will be recorded as the nth run of your algorithm. 
The data that are specified to be recorded for statistical analysis in 
the Technical Report will be recorded automatically in memory.

After all 20 run of the algorithm finished, or part of the running are done, 
you can dump all the recorded result into directory specified as the first 
argument, and all the result fills will have the same prefix as the second 
argument specified. 

The first argument can be in the format like "result/randomsampling", 
"randomsampling", or "d:/result/myaogorithm" under Windows, 
"/home/username/algorithm" under Linux. The main constrains is the 
sub-directory seperator "/" and no "/" for the lowerest level. 
Directory names like "result\\algorithm" or "result/algorithm/" are 
not accepted now.


4. Contact

If you find anything does not work, have any difficulty to use the package, 
or have better ideas about the package, please don't hesitate to contact 
Qin Chen (email:chenqin1980@gmail.com) for further information.
