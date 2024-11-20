function [cons,Data] = Input_Constraints_feedback(Data,V)
%INPUT_CONSTRAINTS generates constraints on the input for every scenario
%realisation V
lower_bound = zeros(Data.Para.Horizon*Data.Para.ul,1);
upper_bound = zeros(Data.Para.Horizon*Data.Para.ul,1);
for i = 0 :Data.Para.Horizon-1
    %Window Blinds
    upper_bound(i*Data.Para.ul+1,1)=1;
    
    % Heating
    upper_bound(i*Data.Para.ul+2,1) = 1e3;
    
    % Cooling
    lower_bound(i*Data.Para.ul+3,1) = -1e3;
end
cons = lower_bound <= Data.Gamma + Data.Theta*V <= upper_bound;
Data.LB = lower_bound;
Data.UB = upper_bound;
end

