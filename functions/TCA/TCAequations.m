function p=TCAequations(q)
global xc yc nx_c ny_c
fai1=q(1);fai2=q(2);beta1=q(3);
p=zeros(3,1);
p(1)=xc*cos(fai2) - yc*sin(fai2) - sin(fai1)*((13*cos(beta1))/2 - 110) + (13*cos(fai1)*sin(beta1))/2;
p(2)=(13*sin(beta1)*sin(fai1))/2 + yc*cos(fai2) + cos(fai1)*((13*cos(beta1))/2 - 110) + xc*sin(fai2) - 5/2;
p(3)=ny_c*cos(fai2) + nx_c*sin(fai2) + (13*cos(beta1)*cos(fai1))/(2*((169*cos(beta1)^2)/4 + (169*sin(beta1)^2)/4)^(1/2)) + (13*sin(beta1)*sin(fai1))/(2*((169*cos(beta1)^2)/4 + (169*sin(beta1)^2)/4)^(1/2));
