function V = Dimcorrect(Data,V_result,k)
%DIMCORRECT this function is used to chagne the dimension of a matrix into
%a vector that is used in the dynamic equation for MPC
V = zeros(Data.Para.Horizon*Data.Para.vl,1);
for i=0 : Data.Para.Horizon-1
    V(i*Data.Para.vl+1:i*Data.Para.vl+Data.Para.vl,1) = V_result(:,i+1,k);
end
end

