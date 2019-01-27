/*
  CEC15 Test Function Suite for Niching Optimization
  Jane Jing Liang 
  email: liangjing@zzu.edu.cn; liangjing@pmail.ntu.edu.cn
  Nov. 25th 2014

  Reference£º	
  B. Y. Qu, J. J. Liang, P. N. Suganthan, Q. Chen, "Problem Definitions and Evaluation Criteria for the CEC 2015 Special Session and Competition on Niching Numerical Optimization",Technical Report201411B,Computational Intelligence Laboratory, Zhengzhou University, Zhengzhou China and Technical Report, Nanyang Technological University, Singapore, November 2014

*/

1. Run "mex -setup" and choose a proper complier (Microsoft Visual C++ 6.0 is prefered)
   Would you like mex to locate installed compilers [y]/n? y
2. Put cec15_nich_func.cpp and input_data folder with your algorithm in the same folder. Set this folder as 	the current path.
3. Run the following command in Matlab window:
   mex cec15_nich_func.cpp -DWINDOWS
4. Then you can use the test functions as the following example:
   f = cec15_nich_func(x,func_num); 
   here x is a D*pop_size matrix.
5. main.m is an example test code with PSO algorithm.
