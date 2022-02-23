function [z pars]=Opt_Monod()
clc

p0=[10      % Initial concentration of glucose
    1e6     % Initial amount of cells
    5e-7    % Vmax for glucose transp
    15      % Km for glucose transp
    1];     % umax max. growth rate

z=Cost(p0,1);

% options=optimset('Display','iter');
% p2 = fminsearch(@Cost,p0,options,0)
% z = Cost(p2,1);
% pars = p2

end

function [z DSE]=Cost(pp,g)
p=abs(pp);

x0=[p(1)
    p(2)];
Tspan=0:.1:10;

[t,y]=ode23s(@My_ode,Tspan,x0,[],p(3:5));

glc=dataset(2,'glc');
cells=dataset(2,'cells');


Sim_glc=interp1(t,y(:,1),glc(:,1));
Sim_cells=interp1(t,y(:,2),cells(:,1));

DSE=[(Sim_glc-glc(:,2)); (Sim_cells-cells(:,2)) ];
Cost=sum(DSE.^2);

if g==1
    subplot(2,1,1)
    plot(t,y(:,2))
    hold on
    plot(cells(:,1),cells(:,2),'o')
    hold off
    xlabel('Time');ylabel('Pop. Size'); title('Simulation of the Monod example')
    subplot(2,1,2)
    plot(t,y(:,1))
    hold on
    plot(glc(:,1),glc(:,2),'o')
    hold off
    xlabel('Time');ylabel('Glucose');
end

z = Cost
DSE = DSE
end

function dxdt=My_ode(t,x,p)

dxdt=[-(p(1)*x(2)*x(1)/(p(2)+x(1)))
    p(3)*x(2)*x(1)/(p(2)+x(1))];

end

function z=dataset(dataset,Met)

switch dataset
    case 2
        switch Met
            case 'glc'
                z=[ 0   40
                    1   39.6
                    2   38.5
                    3   35.4
                    4   27.4
                    5   11.3
                    6    0.7
                    7    0.01
                    8    0
                    9    0
                    10   0];
            case 'cells'
                z=[ 0.25    0.26e6
                    1.25    0.74e6
                    2.25    2.12e6
                    3.25    5.92e6
                    3.75    9.71e6
                    4.25   15.5e6
                    4.75   23.5e6
                    5.25   32.0e6
                    6.25   38.6e6
                    7.25   38.9e6
                   10.00   38.9e6];
        end
        
        
end


end