clear all;
close all;
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

set(groot,'defaultLineLineWidth',1.0);
set(groot,'defaultAxesLineWidth',0.6);

FigH = figure('Position', get(0, 'Screensize'))
hold on
plot(w0,w1_1,'k')
plot(w0,w1_2,'k')
plot(w0_1,w1,'k')
plot(w0_2,w1,'k')
plot(w0_3,w1,'k')
axis([-0.2 0.5 -0.5 0.5])

    set(gca,'xtick',[])
    text(0.48,-0.53,'w_0','Color','k', 'FontSize', 36)
    text(-T0c-0.01,-0.53,'w_0=-T_0^c','Color','b', 'FontSize', 16)
    text(gamma*p*(1/gamma-1)-T0c-0.04,-0.53,{'w_0=\gammap(1/\gamma-1)-T_0^c'},...
         'Color','b', 'FontSize', 16)
    text(gamma*p*(hstar-1)-T0c-0.04,-0.53,'w_0=\gammap(h*-1)-T_0^c',...
         'Color','b', 'FontSize', 16)

    set(gca,'ytick',[])
    text(-0.24,0.48,'w_1','Color','k', 'FontSize', 36)
    text(-0.28,gamma*p*(1-1/gamma)-T1c,'w_1=\gammap(1-1/\gamma)-T_1^c',...
         'Color','b', 'FontSize', 16)
    text(-0.28,gamma*p*(hstar-1/gamma)-T1c,'w_1=\gammap(h*-1/\gamma)-T_1^c',...
         'Color','b', 'FontSize', 16)

    text(0.5*(-2*T0c+gamma*p*(1/gamma-1)),gamma*p*(1-1/gamma)+T0c-T1c+0.3*(-2*T0c+gamma*p*(1/gamma-1)),...
        'w_1=\gammap(1-1/\gamma)+T_0^c-T_1^c+w_0',...
        'Color','b', 'Rotation', 20, 'FontSize', 16)
    text(gamma*p*(hstar-1)-T0c,-T1c,...
         '    w_1=-T_1^c','Color','b', 'FontSize', 16)
    ww0=0.8*gamma*p*(1/gamma-1)+0.2*gamma*p*(hstar-1)-T0c;
    text(ww0+0.01,gamma*p*(hstar-0.25)-T0c-T1c-ww0,...
        'w_1=\gammap(h*-1)-T_0^c-T_1^c-w_0',...
        'Color','b', 'Rotation', -25, 'FontSize', 16)
    clear ww0;
    
    getcoord = @(vec) length(find(w0 <= min(vec))); 
    x = [w0(1), w0(1), w0(getcoord(w0_2)), ...
              w0(getcoord(w0_2)), w0(getcoord(w0_1))];
    y = [w1_1(1), 1.0, 1.0,  w1_1(getcoord(w0_2)), ...
              w1_1(getcoord(w0_1))];
    x2 = [w0(getcoord(w0_2)), w0(getcoord(w0_2)), 1.0, 1.0, ...
              w0(getcoord(w0_3)), w0(getcoord(w0_3))];
    y2 = [w1_1(getcoord(w0_2)), 1.0, 1.0, w1(1), ...
          w1(1), w1_1(getcoord(w0_3))];
    x3 = [w0(getcoord(w0_1)), w0(getcoord(w0_1)), ...
          w0(getcoord(w0_2)), w0(getcoord(w0_3)), w0(getcoord(w0_3))];
    y3 = [w1(1), w1_1(getcoord(w0_1)), w1_1(getcoord(w0_2)), w1_1(getcoord(w0_3)), w1(1)];
    patch(x, y, 'r', 'FaceAlpha', 0.125)
    patch(x2, y2, 'b', 'FaceAlpha', 0.125)
    patch(x3, y3, 'k', 'FaceAlpha', 0.075)

%left side
str={"(DB, DB)"};
%str={'h_0=1','h_1=1'};
text(-T0c-0.091,gamma*p*(1-1/gamma)-T1c-0.05,str, 'FontSize', 24)

str={"(DB, CB)", "h_0 < h_1 < h^*"};
%str={"(DB, CB)", "h_1=1/\gamma+(w_1+T_1^c)/(\gamma p)"};
%str={'h_0=1','h_1=1/\gamma+(w_1+T_1^c)/(\gamma p)'};
text(-T0c-0.091,gamma*p*(0.5*(1-1/gamma)+0.5*(hstar-1/gamma))-T1c+0.1,...
     str, 'FontSize', 24)

str={"(DB, B)", "h_0 < h_1 = h^*"};
%str={"(DB, B)", "h_1=h*"};
%str={'h_0=1','h_1=h*'};
text(-T0c-0.091,(gamma*p*(hstar-1/gamma))-T1c+0.2,str, 'FontSize', 24)

%middle 1
str={"(CB, DB)", "h_{-1} < h_0=h_1 < h^*"};
%str={"(CB, DB)", "h_0=1+(w_0+T_0^c)/(\gamma p)"};
%str={'h_0=1+(w_0+T_0^c)/(\gamma p)','h_1=1+(w_0+T_0^c)/(\gamma p)'};
text(0.4*gamma*p*(1/gamma-1)-T0c,gamma*p*(1-1/gamma)-T1c,str, 'FontSize', 24)

str={"(DB, CB)", "h_0 < h_1 < h^*"};
%str={"(DB, CB)", "h_1=1/\gamma+(w_1+T_1^c)/(\gamma p)"};
%str={'h_0=1','h_1=1/\gamma+(w_1+T_1^c)/(\gamma p)'};
text(0.25*gamma*p*(1/gamma-1)-T0c,gamma*p*(0.5*(1-1/gamma)+0.5*(hstar-1/gamma))-T1c+0.1,...
     str, 'FontSize', 24)

str={"(DB, B)", "h_0 < h_1 = h^*"};
%str={"(DB, B)", "h_1=h*"};
%str={'h_0=1','h_1=h*'};
text(0.25*gamma*p*(1/gamma-1)-T0c,(gamma*p*(hstar-1/gamma))-T1c+0.2,str, 'FontSize', 24)

%middle 2
str={"(CB, DB)", "h_{-1} < h_0=h_1 < h^*"};
%str={"(CB, DB)", "h_0=1+(w_0+T_0^c)/(\gamma p)"};
%str={'h_0=1+(w_0+T_0^c)/(\gamma p)','h_1=h_0'};
text(gamma*p*(0.75*(1/gamma-1)+0.25*(hstar-1))-T0c,0.5*gamma*p*(1-1/gamma)-T1c,str, 'FontSize', 24)

str={"(CB, CB)", "h_{-1} < h_0 < h_1 < h^*"};
%str={"(CB, CB)", "h_t = h_{t-1} + (w_t + T_t^c)/(\gamma p)"};
%str={'h_0=1+(w_0+T_0^c)/(\gamma p)','h_1=h_0+(w_1+T_1^c)/(\gamma p)'};
text(gamma*p*(0.85*(1/gamma-1)+0.10*(hstar-1))-T0c,-T1c+0.375*(gamma*p*(hstar-1)-T0c-(gamma*p*(0.80*(1/gamma-1)+0.20*(hstar-1))-T0c)),...
     str, 'FontSize', 24)

str={"(CB, B)", "h_{-1} < h_0 < h_1=h^*"};
%str={"(CB, B)", "h_0 = 1+(w_0+T_0^c)/(\gamma p)", "h_1=h*"};
%str={'h_0=1+(w_0+T_0^c)/(\gamma p)','h_1=h*'};
text(gamma*p*(0.9*(1/gamma-1)+0.10*(hstar-1))-T0c,(gamma*p*(hstar-1/gamma))-T1c+0.2,...
     str, 'FontSize', 24)

%right
str={"(B, B)", "h_0=h_1=h^*"};
%str={'h_0=h*','h_1=h*'};
text(gamma*p*(hstar-1)-T0c+0.05,(gamma*p*(hstar-1/gamma))-T1c+0.2,...
     str, 'FontSize', 24)

h = zeros(3, 1);
h(1) = patch(NaN,NaN,'sk', 'FaceAlpha', 0.15);
h(2) = patch(NaN,NaN,'sr', 'FaceAlpha', 0.15);
h(3) = patch(NaN,NaN,'sb', 'FaceAlpha', 0.15);
legend(h, {'Purchase in Period 0','Purchase in Period 1','Purchase in Both Periods'},...
       'Location', 'southeast', 'FontSize', 28);

    set(gcf,'color','w');
    %suptitle('Housing Choice Thresholds')
    savefig('Thresholds')
    F    = getframe(FigH);
    imwrite(F.cdata,'Thresholds.png', 'png')
    %close all;
    %clear fig F FigH;



T0c2=0.15;
w0_1(2,:)=-T0c2*ones(size(w1));
w0_2(2,:)=(gamma*p*(1/gamma-1)-T0c2)*ones(size(w1));
w0_3(2,:)=(gamma*p*(hstar-1)-T0c2)*ones(size(w1));

w1_1(2,:)=-T1c+(gamma*p*(1-1/gamma))*(w0<=-T0c2)+(gamma*p*(1-1/gamma)+T0c2+w0).*(w0>-T0c2).*(w0<=(gamma*p*(1/gamma-1)-T0c2));
w1_1(2,w0>gamma*p*(hstar-1)-T0c2)=nan;

w1_2(2,:)=(gamma*p*(hstar-1/gamma)-T1c)*(w0<=(gamma*p*(1/gamma-1)-T0c2))+(gamma*p*(hstar-1)-T1c-T0c2-w0).*(w0>(gamma*p*(1/gamma-1)-T0c2));
w1_2(2,w0>gamma*p*(hstar-1)-T0c2)=nan;

FigH = figure('Position', get(0, 'Screensize'))  
hold on
plot(w0,w1_1(1,:),'--','Color','k')
plot(w0,w1_2(1,:),'--','Color','k')
plot(w0_1(1,:),w1,'--','Color','k')
plot(w0_2(1,:),w1,'--','Color','k')
plot(w0_3(1,:),w1,'--','Color','k')
plot(w0,w1_1(2,:),'k')
plot(w0,w1_2(2,:),'k')
plot(w0_1(2,:),w1,'k')
plot(w0_2(2,:),'k')
plot(w0_3(2,:),w1,'k')
axis([-0.2 0.5 -0.5 0.5])
text(0.4,0.4,{'Dashed line: T_0^c=0.1','Solid line: T_0^c=0.15'})
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gcf,'color','w');
    %suptitle('Housing Choice Thresholds Change to T_0^c\uparrow')
    savefig('ThresholdsChangeT0c')
    F    = getframe(FigH);
    imwrite(F.cdata,'ThresholdsChangeT0c.png', 'png')
    %close all;
    %clear fig F FigH;

