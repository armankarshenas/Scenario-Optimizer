function [Data] = Season(Data,Summer)
if Summer ==1
    Data.Para.Tmax = 24;
    Data.xi = 20*ones(Data.Para.xl,1);
    Data.Para.Amb = 35;
    Data.Para.Tmin =0;
elseif Summer ==0
    Data.Para.Tmin = 15;
    Data.xi = 18*ones(Data.Para.xl,1);
    Data.Para.Amb = 5;
    Data.Para.Tmax =0;
else
    fprintf('Season must be selected');
end
end

