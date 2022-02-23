function z=Simple_Opt()
clc

p0=[1   % the lower asymptote
    1   % the upper asymptote
    1   % the growth rate
    1   % how max. growth occurs
    1   % Depend. on y(0)
    1]; % the time of max growth

z=Cost(p0,1);

% p2=fminsearch(@Cost,p0,[],0)
% z=Cost(p2,1);


end

function z=Cost(p,g)

Exp_data=dataset(1);

t=[0:.1:10]';
y=p(1)+(p(2)-p(1))./(1+p(5)*exp(-p(3)*(t-p(6)))).^(1/p(4));

ys=interp1(t,y,Exp_data(:,1));

Cost=sum((Exp_data(:,2)-ys).^2);

if g==1
    subplot(2,1,1)
    plot(t,y)
    hold on
    plot(Exp_data(:,1),Exp_data(:,2),'o')
    hold off
    xlabel('Time');ylabel('Pop Size'); title('Simulation of the generalised logistic function')
    subplot(2,1,2)
    plot(Exp_data(:,1),Exp_data(:,2)-ys,'o')
    xlabel('Time');ylabel('Residuals'); title('Residuals Plot')
end

z=Cost;
end

function z=dataset(dataset)

switch dataset
    case 1
        z=[      0   12.0021
            0.5000   12.0047
            1.0000   12.0109
            1.5000   12.0251
            2.0000   12.0577
            2.5000   12.1328
            3.0000   12.3056
            3.5000   12.7032
            4.0000   13.6180
            4.1000   13.9114
            4.4000   15.1509
            4.7000   17.1909
            5.0000   20.5287
            5.3000   25.8510
            5.6000   33.5000
            5.9000   41.4622
            6.0000   43.4489
            6.5000   47.5205
            7.0000   47.9597
            7.5000   47.9967
            8.0000   47.9997
            8.5000   48.0000
            9.0000   48.0000
            9.5000   48.0000
            10.0000   48.0000];
end


end