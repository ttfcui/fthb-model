clear all;
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

w1_3(1,:)=(gamma*p*(hstar-1)-T1c-T0c-w0).*(w0>=-T0c);
w1_3(1,w0<-T0c)=nan;
w1_3(1,w0>gamma*p*(hstar-1)-T0c)=nan;
w1_3(1,:) = nan;


set(groot,'defaultLineLineWidth',1.0);
set(groot,'defaultAxesLineWidth',0.6);

fig= figure;
plot(w0,w1_1,'k',w0,w1_2,'k',w0,w1_3,'k',w0_1,w1,'k',w0_2,w1,'k',w0_3,w1,'k')
%legend('w_{1,1}','w_{1,2}','w_{1,3}','w_{0,1}','w_{0,2}','w_{0,3}','Location','Northwest')
axis([-0.2 0.5 -0.5 0.5])
%close all;


T0c2=0.15;
T1c = .10;
gamma = .10;
w0_1(2,:)=-T0c2*ones(size(w1));
w0_2(2,:)=(gamma*p*(1/gamma-1)-T0c2)*ones(size(w1));
w0_3(2,:)=(gamma*p*(hstar-1)-T0c2)*ones(size(w1));

w1_1(2,:)=-T1c+(gamma*p*(1-1/gamma))*(w0<=-T0c2)+(gamma*p*(1-1/gamma)+T0c2+w0).*(w0>-T0c2).*(w0<=(gamma*p*(1/gamma-1)-T0c2));
w1_1(2,w0>gamma*p*(hstar-1)-T0c2)=nan;

w1_2(2,:)=(gamma*p*(hstar-1/gamma)-T1c)*(w0<=(gamma*p*(1/gamma-1)-T0c2))+(gamma*p*(hstar-1)-T1c-T0c2-w0).*(w0>(gamma*p*(1/gamma-1)-T0c2));
w1_2(2,w0>gamma*p*(hstar-1)-T0c2)=nan;

w1_3(2,:)=(gamma*p*(hstar-1)-T1c-T0c2-w0).*(w0>=-T0c2);
w1_3(2,w0<-T0c2)=nan;
w1_3(2,w0>gamma*p*(hstar-1)-T0c2)=nan;
w1_3(2,:) = nan;

plot(w0,w1_1(1,:),'--k',w0,w1_1(2,:),'-k',...
     w0,w1_2(1,:),'--k',w0,w1_2(2,:),'-k',...
     w0,w1_3(1,:),'--k',w0,w1_3(2,:),'-k',...
     w0_1(1,:),w1,'--k',w0_1(2,:),w1,'-k',...
     w0_2(1,:),w1,'--k',w0_2(2,:),w1,'-k',...
     w0_3(1,:),w1,'--k',w0_3(2,:),w1,'-k')
%legend('w_{1,1}','w_{1,2}','w_{1,3}','w_{0,1}','w_{0,2}','w_{0,3}','Location','Northwest')
line(w0, w0-(T0c2-T1c), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2.5)
    text(-0.075,-0.155,...
        'f(w_0) = w_0+ (T_1^c - T_0^c)',...
        'Color','b', 'Rotation', 25, 'FontSize', 12)

axis([-0.2 0.5 -0.5 0.5])

getcoord = @(vec) length(find(w0 <= min(vec)));
% Pos extensive, (DB, DB) -> (CB, DB)
x_ext1 = [w0(getcoord(w0_1(2,:))), w0(getcoord(w0_1(2,:))), ...
          w0(getcoord(w0_1(1,:))), w0(getcoord(w0_1(1,:)))];
y_ext1 = [w1(1), w1_1(1,getcoord(w0_1(1,:))),...
          w1_1(1,getcoord(w0_1(1,:))), w1(1)];
% Pos extensive, (DB, CB/B) -> (CB, CB/B)
x_ext2 = [w0(getcoord(w0_2(2,:))), w0(getcoord(w0_2(2,:))), ...
          w0(getcoord(w0_2(1,:))), w0(getcoord(w0_2(1,:)))];
y_ext2 = [w1_1(2, getcoord(w0_2(2,:))), 1.0, 1.0, w1_1(2, getcoord(w0_3(2,:)))];
% Pos extensive, (CB, DB) -> (CB, CB/B)
x_ext3 = [w0(getcoord(w0_2(1,:))), w0(getcoord(w0_2(1,:))), ...
          w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(2,:)))];
y_ext3 = [w1_1(2, getcoord(w0_2(1,:))), ...
          max(w1_1(2, getcoord(w0_2(1,:))), w1_1(1, getcoord(w0_2(1,:)))),...
          max(w1_1(2, getcoord(w0_2(1,:))), w1_1(1, getcoord(w0_2(1,:)))),...
          w1_1(2, getcoord(w0_3(2,:)))];
% Pos extensive, (B, DB) -> (CB, CB/B)
x_ext4 = [w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(2,:))), ...
      min(w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(1,:)))),...
      min(w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(1,:))))];
y_ext4 = [w1_2(1, getcoord(w0_3(1,:))), 1.0, 1.0, w1_2(1, getcoord(w0_3(1,:)))];
% Pos extensive, (DB, DB) -> (DB, CB)
x_ext5 = [w0(1), w0(getcoord(w0_1(2,:))), w0(getcoord(w0_1(2,:))), w0(1)];
y_ext5 = [w0(getcoord(w1_1(1,:))), w0(getcoord(w1_1(1,:))),...
          min(w0(getcoord(w1_1(1,:))), w0(getcoord(w1_1(2,:)))),...
          min(w0(getcoord(w1_1(1,:))), w0(getcoord(w1_1(2,:))))];
% Pos extensive, (CB, DB) -> (B, B)
x_ext6 = [w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(2,:))), ...
      w0(getcoord(w0_3(1,:))), w0(getcoord(w0_3(1,:)))];
y_ext6 = [w1(1), w1_1(1, getcoord(w0_3(2,:))), ...
          min(w1_1(1, getcoord(w0_3(1,:))), w1_1(2, getcoord(w0_3(1,:)))),...
          w1(1)];
patch(x_ext1, y_ext1, 'r', 'FaceAlpha', 0.125)
patch(x_ext2, y_ext2, 'r', 'FaceAlpha', 0.125)
patch(x_ext3, y_ext3, 'r', 'FaceAlpha', 0.125)
patch(x_ext4, y_ext4, 'r', 'FaceAlpha', 0.125)
patch(x_ext5, y_ext5, 'r', 'FaceAlpha', 0.125)
patch(x_ext6, y_ext6, 'r', 'FaceAlpha', 0.125)

% Inframarginal upscaling, (CB, CB) -> (CB, B)
area_intersect = length(find(w1_2(2,:) >= w1_2(1, getcoord(w0_3(1,:)))));
x_int1 = [w0(getcoord(w0_2(1,:))), w0(getcoord(w0_2(1,:))), ...
          w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(2,:))), w0(area_intersect)];
y_int1 = [w1_2(2, getcoord(w0_2(1,:))), w1_2(1, getcoord(w0_2(1,:))),...
          w1_2(1, getcoord(w0_3(2,:))), w1_2(2, area_intersect),...
          w1_2(2, area_intersect)];
% Inframarginal upscaling, (CB, DB) -> (B, DB)
x_int2 = [w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(2,:))), ...
      max(w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(1,:)))),...
      max(w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(1,:))))];
y_int2 = [w1_2(1, getcoord(w0_3(1,:))), 1.0, 1.0, w1_2(1, getcoord(w0_3(1,:)))];
% Inframarginal upscaling, (DB, CB) -> (DB, B)
x_int3 = [w0(1), w0(getcoord(w0_2(2,:))), w0(getcoord(w0_2(2,:))), ...
          w0(1)];
y_int3 = [w1_2(2, 1), w1_2(2, getcoord(w0_2(2,:))),...
          max(w1_2(1, getcoord(w0_2(2,:))), w1_2(2, getcoord(w0_2(2,:)))),...
          max(w1_2(1, getcoord(w0_2(2,:))), w1_2(2, getcoord(w0_2(2,:))))];
patch(x_int1, y_int1, 'k', 'FaceAlpha', 0.075)
patch(x_int2, y_int2, 'k', 'FaceAlpha', 0.075)
patch(x_int3, y_int3, 'k', 'FaceAlpha', 0.075)

% Neg extensive, (B, B) -> (CB, DB)
x3 = [w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(2,:))), ...
      min(w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(1,:)))),...
      min(w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(1,:))))];
y3 = [w1(1), w1_1(1, getcoord(w0_3(2,:))), ...
      max(w1_1(1, getcoord(w0_3(1,:))), w1_1(2, getcoord(w0_3(1,:)))),...
          w1(1)];
% Neg extensive, (CB, CB/B) -> (CB, DB)
x32 = [w0(getcoord(w0_2(1,:))), w0(getcoord(w0_2(1,:))), ...
          w0(getcoord(w0_3(2,:))), w0(getcoord(w0_3(2,:)))];
y32 = [w1_1(2, getcoord(w0_2(1,:))), ...
       min(w1_1(2, getcoord(w0_2(1,:))), w1_1(1, getcoord(w0_2(1,:)))),...
          min(w1_1(2, getcoord(w0_2(1,:))), w1_1(1, getcoord(w0_2(1,:)))),...
          w1_1(2, getcoord(w0_3(2,:)))];
% Neg extensive, (DB, CB) -> (DB, DB)
x33 = [w0(1), w0(getcoord(w0_1(2,:))), w0(getcoord(w0_1(2,:))), w0(1)];
y33 = [w0(getcoord(w1_1(1,:))), w0(getcoord(w1_1(1,:))),...
          max(w0(getcoord(w1_1(1,:))), w0(getcoord(w1_1(2,:)))),...
          max(w0(getcoord(w1_1(1,:))), w0(getcoord(w1_1(2,:))))];
patch(x3, y3, 'y', 'FaceAlpha', 0.15)
patch(x32, y32, 'y', 'FaceAlpha', 0.15)
patch(x33, y33, 'y', 'FaceAlpha', 0.15)


x4 = [w0(getcoord(w0_1(2,:))), w0(getcoord(w0_2(2,:))), ...
      w0(getcoord(w0_2(1,:))), w0(getcoord(w0_2(1,:))), ...
      w0(getcoord(w0_1(1,:))), w0(getcoord(w0_1(2,:)))];
y4 = [w1_1(2, getcoord(w0_1(2,:))), w1_1(2, getcoord(w0_2(2,:))),...
      max(w1_1(1, getcoord(w0_2(1,:))), w1_1(2, getcoord(w0_2(1,:)))),...
      w1_1(1, getcoord(w0_2(1,:))), w1_1(1, getcoord(w0_1(1,:))), ...
      min(w1_1(1, getcoord(w0_1(2,:))), w1_1(2, getcoord(w0_1(2,:))))];
patch(x4, y4, 'b', 'FaceAlpha', 0.15)

set(gcf,'color','w');
suptitle('T_0^c = 0.15, T_1^c = 0.10, \gamma=0.1');


figure;
h = zeros(4, 1);
h(1) = patch(NaN,NaN,'sk', 'FaceAlpha', 0.15);
h(2) = patch(NaN,NaN,'sr', 'FaceAlpha', 0.15);
h(3) = patch(NaN,NaN,'sy', 'FaceAlpha', 0.15);
h(4) = patch(NaN,NaN,'sb', 'FaceAlpha', 0.15);
legend(h, {'Intensive Margin','Positive Extensive Margin','Negative Extensive Margin'...
       'Timing Margin'}, 'Location', 'eastoutside', 'FontSize', 16,...
       'NumColumns', 2);
set(gcf,'color','w');