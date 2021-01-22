w0=-1:0.00001:1;
w1=w0;

gamma=0.1;
p=1/3;
hstar=15;

T0c=0.1;
T1c=0.1;

%% If T0b does not matter, set it to zero, this is what we get

w0_1(1,:)=-T0c*ones(size(w1));
w0_2(1,:)=(gamma*p*(1/gamma-1)-T0c)*ones(size(w1));
w0_3(1,:)=(gamma*p*(hstar-1)-T0c)*ones(size(w1));

w1_1(1,:)=-T1c+(gamma*p*(1-1/gamma))*(w0<=-T0c)+(gamma*p*(1-1/gamma)+T0c+w0).*(w0>-T0c).*(w0<=(gamma*p*(1/gamma-1)-T0c));
w1_1(1,w0>gamma*p*(hstar-1)-T0c)=nan;

w1_2(1,:)=(gamma*p*(hstar-1/gamma)-T1c)*(w0<=(gamma*p*(1/gamma-1)-T0c))+(gamma*p*(hstar-1)-T1c-T0c-w0).*(w0>(gamma*p*(1/gamma-1)-T0c));
w1_2(1,w0>gamma*p*(hstar-1)-T0c)=nan;
%w1_2(1,w0<=(gamma*p*(1/gamma-1)-T0c))=nan;

%w1_3(1,:)=(gamma*p*(hstar-1)-T1c-T0c-w0).*(w0>=-T0c);
w1_3(1,w0<-T0c)=nan;
w1_3(1,w0>gamma*p*(hstar-1)-T0c)=nan;
w1_3(1,:) = nan;


set(groot,'defaultLineLineWidth',1.0);
set(groot,'defaultAxesLineWidth',0.6);

plot(w0,w1_1,'k',w0,w1_2,'k',w0,w1_3,'k',w0_1,w1,'k',w0_2,w1,'k',w0_3,w1,'k')
%legend('w_{1,1}','w_{1,2}','w_{1,3}','w_{0,1}','w_{0,2}','w_{0,3}','Location','Northwest')
axis([-0.2 0.5 -0.5 0.5])
%close all;
