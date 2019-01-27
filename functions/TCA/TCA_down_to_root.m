function X1 =TCA_down_to_root( initial_point,DPoint,U,U_d,Q1,BasicParameter1 )
% TCA_down_to_root.m  用牛顿迭代法求解TCA方程，往齿根计算
% initial_point       输入初始参考点的u/beta/fai10,fai20
% DPointU,U_d,Q1      输入拟合基本参量,分别为控制顶点，节点矢量，导函数节点矢量，导函数控制顶点
% BasicParameter1     输入摆线针轮基本参量
% X1                  输出TCA方程的解


tol=0.000001;
%摆阵针轮基本参数
zc=BasicParameter1(1);
zp=BasicParameter1(2);
rp=BasicParameter1(3);
a=BasicParameter1(4);
rrp=BasicParameter1(5);
d_rrp=BasicParameter1(6);
d_rp=BasicParameter1(7);
deta=BasicParameter1(8);




%初始参考点
u0=initial_point(1);
beta10=initial_point(2);
fai1_st=initial_point(3);
fai20=initial_point(4);
p=3;


m1=1;
X1=zeros(2,4);
for fai1=fai1_st:0.01:3.142
    N=100;
    fai10=fai1;
    for i=1:N
        [ f_rc1,f_nc1 ] = Surfacepoint_BYT( u0,p,U,DPoint,U_d,Q1,fai20,a);
        [ f_rp1,f_np1 ] = Surfacepoint_Pin(BasicParameter1,beta10,fai10);
               
        F1=f_rc1(1)-f_rp1(1);
        F2=f_rc1(2)-f_rp1(2);
        F3=f_nc1(2)-f_np1(2);
        
        FN=[F1,F2,F3];
        F_distance=norm(FN);
        if F_distance < tol
            u_output=u0;
            beta1_output=beta10;
            fai1_output=fai10;
            fai2_output=fai20;
           
            X1 (m1,:)=[ u_output,beta1_output,fai1_output,fai2_output];
            break;
        end
            
     B=[-F1,-F2,-F3]';
     dj=0.000001;
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     u_plus=u0+dj;
     [ f_rc11,f_nc11 ] = Surfacepoint_BYT( u_plus,p,U,DPoint,U_d,Q1,fai20,a);
     [ f_rp11,f_np11 ] = Surfacepoint_Pin(BasicParameter1,beta10,fai10);
      Fu11=f_rc11(1)-f_rp11(1);
      Fu21=f_rc11(2)-f_rp11(2);
      Fu31=f_nc11(2)-f_np11(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      u_sub=u0-dj;
      [ f_rc12,f_nc12 ] = Surfacepoint_BYT( u_sub,p,U,DPoint,U_d,Q1,fai20,a);
      [ f_rp12,f_np12 ] = Surfacepoint_Pin(BasicParameter1,beta10,fai10);
      Fu12=f_rc12(1)-f_rp12(1);
      Fu22=f_rc12(2)-f_rp12(2);
      Fu32=f_nc12(2)-f_np12(2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
      dFu1=(Fu11-Fu12)/(2*dj);
      dFu2=(Fu21-Fu22)/(2*dj);
      dFu3=(Fu31-Fu32)/(2*dj);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
       beta1_plus=beta10+dj;
       [ f_rc21,f_nc21 ] = Surfacepoint_BYT( u0,p,U,DPoint,U_d,Q1,fai20,a);
       [ f_rp21,f_np21 ] = Surfacepoint_Pin(BasicParameter1,beta1_plus,fai10);
       Fbet11=f_rc21(1)-f_rp21(1);
       Fbet21=f_rc21(2)-f_rp21(2);
       Fbet31=f_nc21(2)-f_np21(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
       beta1_sub=beta10-dj;
       [ f_rc22,f_nc22 ] = Surfacepoint_BYT( u0,p,U,DPoint,U_d,Q1,fai20,a);
       [ f_rp22,f_np22 ] = Surfacepoint_Pin(BasicParameter1,beta1_sub,fai10);
       Fbet12=f_rc22(1)-f_rp22(1);
       Fbet22=f_rc22(2)-f_rp22(2);
       Fbet32=f_nc22(2)-f_np22(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
      dFbet1=(Fbet11-Fbet12)/(2*dj);
      dFbet2=(Fbet21-Fbet22)/(2*dj);
      dFbet3=(Fbet31-Fbet32)/(2*dj); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     fai2_plus=fai20+dj;
     [ f_rc31,f_nc31 ] = Surfacepoint_BYT( u0,p,U,DPoint,U_d,Q1,fai2_plus,a);
     [ f_rp31,f_np31 ] = Surfacepoint_Pin(BasicParameter1,beta10,fai10);
      Ffai2_11=f_rc31(1)-f_rp31(1);
      Ffai2_21=f_rc31(2)-f_rp31(2);
      Ffai2_31=f_nc31(2)-f_np31(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      fai2_sub=fai20-dj;
      [ f_rc32,f_nc32 ] = Surfacepoint_BYT( u0,p,U,DPoint,U_d,Q1,fai2_sub,a);
      [ f_rp32,f_np32 ] = Surfacepoint_Pin(BasicParameter1,beta10,fai10);
      Ffai2_12=f_rc32(1)-f_rp32(1);
      Ffai2_22=f_rc32(2)-f_rp32(2);
      Ffai2_32=f_nc32(2)-f_np32(2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
      dFfai2_1=(Ffai2_11-Ffai2_12)/(2*dj);
      dFfai2_2=(Ffai2_21-Ffai2_22)/(2*dj);
      dFfai2_3=(Ffai2_31-Ffai2_32)/(2*dj);   
%       Jacobi_ubetfai2=[dFu1,dFu2,dFu3;
%                   dFbet1,dFbet2,dFbet3;
%                   dFfai2_1,dFfai2_2,dFfai2_3];
      Jacobi_ubetfai2=[dFu1,dFbet1,dFfai2_1;
                       dFu2,dFbet2,dFfai2_2;
                       dFu3,dFbet3,dFfai2_3;];
     xk= Jacobi_ubetfai2\B; 
      uq=u0+xk(1);
      beta1q=beta10+xk(2);
      fai2q=fai20+xk(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      if uq < 0 || uq>1
         disp('已在齿根退出啮合') 
        return;
      end
      
      q0=[u0,beta10,fai20];
      q=[uq,beta1q,fai2q];
      u0=q(1);
      beta10=q(2);
      fai20=q(3);
      
      err=norm(q-q0);
      if (err<tol)
        break;
      end
    end
   X1(m1,:)=[u0,beta10,fai10,fai20];
    m1=m1+1;
end
end

