function [V_forecast] = forecast(Data)
% This function generates a reference vector of uncertainty used for
% comparison in all three approaches
V_forecast = zeros(Data.Para.Horizon*Data.Para.vl,1);
P_rand = Data.Para.Lambda;
for i=0:Data.Para.Horizon-1
    change = poissrnd(0.1*Data.Para.Lambda,1);
    % Applying a compact support
    if change >2
        change =2;
    end
    red_inc = rand;
    coef = abs(rand)*(0.948-0.606)+0.606;
    V_forecast(i*Data.Para.vl+1:i*Data.Para.vl+Data.Para.vl,1) = [(P_rand+sign(red_inc)*change)*coef*12; Data.Para.Amb; Data.Para.GR];
end
end

