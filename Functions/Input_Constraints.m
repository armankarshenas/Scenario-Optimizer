function [Constraints] = Input_Constraints(Data)
%INPUT_CONSTRAINTS generates the constraints on the input
Constraints = [];
% Constraints on the blinds
for i = 0 :Data.Para.Horizon-1
    Constraints = [Constraints, 0<=Data.Gamma(i*Data.Para.ul+1)<=1];
end
% Constraints on the heater
for i=0 : Data.Para.Horizon-1
    Constraints = [Constraints ,  0<=Data.Gamma(i*Data.Para.ul+2)<= 1e3];
end
% Constraints on the cooler

for i=0 : Data.Para.Horizon-1
    Constraints =  [Constraints , -1e3<= Data.Gamma(i*Data.Para.ul+3)<=0];
end
end

