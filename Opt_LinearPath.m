function [z pars]=Opt_LinearPath()
clc

p0=[1     % Rate constant for V1
    0.5   % Kinetic order for V1
    1     % Rate constant for V2
    .5    % Kinetic order for V2
    1     % Rate constant for V3
    .5];  % Kinetic order for V3

z=Cost(p0,1);

% options=optimset('Display','iter');
% p2=fminsearch(@Cost,p0,options,0)
% z=Cost(p2,1);
% pars = p2

end

function z=Cost(pp,g)
p=abs(pp);

x0=[10
    .001
    .001
    .001];
Tspan=0:.05:8;

[t,y]=ode23s(@My_ode,Tspan,x0,[],p);

DataE=dataset(3);

DataS=interp1(t,y(:,:),DataE(:,1));

DSE=[DataE(:,2:5)-DataS];

Cost=sum(sum(DSE.^2));

if g==1
    subplot(2,2,1)
    plot(t,y(:,1))
    hold on
    plot(DataE(:,1),DataE(:,2),'+')
    hold off
    xlabel('Time');ylabel('Substrate (X1)'); title('Extracellular')
    subplot(2,2,2)
    plot(t,y(:,4))
    hold on
    plot(DataE(:,1),DataE(:,5),'+')
    hold off
    xlabel('Time');ylabel('Product (X4)'); title('Extracellular')
    subplot(2,2,3)
    plot(t,y(:,2))
    hold on
    plot(DataE(:,1),DataE(:,3),'*')
    hold off
    xlabel('Time');ylabel('Intermediate (X2)');  title('Intracellular')
    subplot(2,2,4)
    plot(t,y(:,3))
    hold on
    plot(DataE(:,1),DataE(:,4),'*')
    hold off
    xlabel('Time');ylabel('Intermediates (X3)');  title('Intracellular')
end

z=Cost;
end

function dxdt=My_ode(t,x,p)

Vi=2.5; Ve=50;


 v=[p(1)*x(1)^p(2)
    abs(p(3))*x(2)^abs(p(4))
    abs(p(5))*x(3)^abs(p(6))];

dxdt=[-v(1)/Ve
       v(1)/Vi-v(2)/Vi
       v(2)/Vi-v(3)/Vi
       v(3)/Ve];

end

function z=dataset(dataset,Met)

switch dataset

        case 3
        z=[0   10.0000    0.0010    0.0010    0.0010
    0.1000    9.4838    5.3635    1.5358    0.1723
    0.2000    9.0006    7.0424    3.6830    0.4643
    0.3000    8.5477    7.3807    5.4507    0.8119
    0.4000    8.1229    7.1486    6.6203    1.1898
    0.5000    7.7242    6.6813    7.1927    1.5832
    0.6000    7.3496    6.1415    7.2537    1.9817
    0.7000    6.9975    5.6013    6.9140    2.3779
    0.8000    6.6662    5.0926    6.2931    2.7656
    0.9000    6.3543    4.6263    5.5022    3.1404
    1.0000    6.0604    4.2041    4.6388    3.4986
    1.2000    5.5218    3.4830    3.0099    4.1546
    1.4000    5.0418    2.8999    1.8135    4.7236
    1.6000    4.6129    2.4259    1.1013    5.2119
    1.8000    4.2286    2.0386    0.7046    5.6354
    2.0000    3.8835    1.7205    0.4693    6.0081
    2.5000    3.1630    1.1454    0.1854    6.7716
    3.0000    2.6021    0.7795    0.0789    7.3561
    3.5000    2.1602    0.5410    0.0354    7.8121
    4.0000    1.8082    0.3820    0.0166    8.1729
    5.0000    1.2953    0.1994    0.0040    8.6956
    6.0000    0.9521    0.1097    0.0011    9.0435
    7.0000    0.7153    0.0630    0.0003    9.2826
    8.0000    0.5478    0.0376    0.0001    9.4514];
        
end


end