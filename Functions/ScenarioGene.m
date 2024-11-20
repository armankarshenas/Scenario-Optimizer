function [V,v_local] = ScenarioGene(Data)
%SCENARIOGENE produces vector of dimension v_length*M that represent a
%scenario Vector V is the v_lengthM * 1 vector used for inequality
%constraints and v_local is the v_length* M matrix
v_length = Data.Para.vl;
v_local = zeros(3,Data.Para.Horizon);
P_rand = poissrnd(Data.Para.Lambda,1);
% Compact Support
if P_rand > 5
    P_rand = 5;
end
for i=1:Data.Para.Horizon
    change = poissrnd(0.1*Data.Para.Lambda,1);
    % Compact Support
    if change >2
        change =2;
    end
    red_inc = rand;
    coef = abs(rand)*(0.948-0.606)+0.606;
    v_local(:,i) = [(P_rand+sign(red_inc)*change)*coef*12 Data.Para.Amb Data.Para.GR];
end
for i=0:(Data.Para.Horizon-1)
    V(i*v_length+1:i*v_length+v_length,1) = v_local(:,(i+1));
end

end

