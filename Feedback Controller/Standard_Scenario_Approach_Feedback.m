%% Scenario approach in MPC for energy management in a building
%% Initation and Dynamic constraint
close all;
clear;
clc;
addpath(genpath('/Users/Karshenas/Arman/Research + Paper/2019-2020/EUROP 2019/tbxmanager/toolboxes'));
Sampling = 0.25;
Data.Para.Sampling = Sampling;

% MPC Horizon and scenarios
Horizon_time = 2;
M = Horizon_time/Sampling;
Data.Para.Horizon = M;
Data.Para.Horizon_Time = Horizon_time;

% Matrix Generation
[F,G,H,Q,R,Z,Data] = MatrixGene(Data);

% Identifiers for x,u and v
x_length = length(Data.Building.building_model.identifiers.x);
v_length = length(Data.Building.building_model.identifiers.v);
u_length = length(Data.Building.building_model.identifiers.u);
Data.Para.xl = x_length;
Data.Para.ul = u_length;
Data.Para.vl = v_length;

% The input variable / decision vector for u_radiator, u_heater, u_blinds
Gamma = sdpvar(M*u_length,1);
Theta = sdpvar(M*u_length,M*v_length);
Data.Gamma = Gamma;
Data.Theta = Theta;
% Summer/Winter:
Summer =1;
Data = Season(Data,Summer);
Data.Summer = Summer;
% The state vector for temperatures
x_o = Data.xi;

% Optimal Temperature Range
t_min = Data.Para.Tmin;
t_max = Data.Para.Tmax;

% The temperature range for comfort
g_t_min = zeros(M*x_length,1);
g_t_max = 100*ones(M*x_length,1);
for i=0:M-1
    g_t_min(i*x_length+1:i*x_length+3) = [t_min ; t_min; t_min];
    g_t_max(i*x_length+1:i*x_length+3) = [t_max ; t_max; t_max];
end

Data.Tmin = g_t_min;
Data.Tmax = g_t_max;
% Radiation and Convection range by a person
h_min = 0.844;
h_max = 0.948;

% Scenario Approach Initiation
epsilon = 0.1;
beta = 1e-4;
d = M*u_length;
N = ceil((2/epsilon)*(log(1/beta)+d));
Data.Para.Eps = epsilon;
Data.Para.Beta = beta;

% Matrices to store the results
X_result = zeros(x_length,M);
V_result = zeros(v_length,M);
U_result = zeros(u_length,M);

% Testing Scenarios parametes
Ntest = 3000;
N_run = 100;
% The uncertainty vector parameters
lambda  = 3;
Amb_temp = Data.Para.Amb;
Global_radiation = 200;
Data.Para.Lambda = lambda;
Data.Para.GR = Global_radiation;

expected =Data.Para.Lambda*ones(M*v_length,1)*9.324;
secondmoment = (Data.Para.Lambda*9.324+(Data.Para.Lambda*9.324)^2)*ones(M*v_length,M*v_length);

%% Constraints on Theta
Constraints = [];
zero_theta = zeros(u_length,v_length);
for i=0:M-1
    for j=i:M-1
        Constraints = [ Constraints Theta(i*u_length+1:i*u_length + u_length,j*v_length+1:j*v_length+v_length) == zero_theta];
    end
end
%% Scenarios Generating Constraints
maxv=zeros(M,1);

for k=1:N
    [V,v_local] = ScenarioGene(Data);
    V_result(:,:,k) =v_local;
    if Summer ==1
        Constraintsinput = Input_Constraints_feedback(Data,V);
        Constraints =  [Constraints , F*x_o+G*Gamma+(H+G*Theta)*V <= g_t_max , Constraintsinput;];
        
    elseif Summer ==0
        Constraintsinput = Input_Constraints_feedback(Data,V);
        Constraints = [Constraints , F*x_o+G*Gamma+(H+G*Theta)*V >= g_t_min, Constraintsinput;];
    else
        fprintf('Season must be selected');
    end
    
end
%{
minv=V(M*v_length-2,1)*ones(M,1);
for i=1:M
    for k=1:N
        V = Dimcorrect(Data,V_result,k);
        if V((i-1)*u_length+1,1) > maxv(i)
            maxv(i) = V((i-1)*u_length+1,1);
        end
        if V((i-1)*u_length+1,1) < minv(i)
            minv(i) = V((i-1)*u_length+1,1);
        end
        
    end
end
Vmax = V;
for i=0:M-1
    Vmax(i*u_length+1,1) = maxv(i+1);
end
Vmin = V;
for i=0:M-1
    Vmin(i*u_length+1,1) = minv(i+1);

%}


%% Cost Function
Cost = Gamma'*R*Gamma + 2*Gamma'*R*Theta*expected+trace(Theta'*R*Theta*secondmoment);
%% Scenarios and Optimisation
options = sdpsettings('solver','quadprog','debug',1);
solution = optimize(Constraints,Cost,options);
U = value(Gamma)+value(Theta)*forecast(Data);
Data.U = U;
for i=0:M-1
    U_result(:,i+1) = U(i*u_length+1:i*u_length+u_length,1);
end
% Generating the states
X = F*x_o+G*U+H*V;
for i=0:M-1
    X_result(:,i+1) = X(i*x_length+1:i*x_length+x_length,1);
end
Data.Xr = X_result;
Data.Ur = U_result;
Data.Vr = V_result;
%% Plotting Initiation
t = 0 : Sampling : Horizon_time;
vec_t_min = t_min*ones(1,M+1);
vec_t_max = t_max*ones(1,M+1);
figure; title('Temperature of different zones / C')
%% Violation Counter for Scenario Approach
violation_counter = violation(Data,Ntest,1)
violation_ratio = violation_counter/Ntest
%% Plotting the results
Plot_Optimal(Data,t,vec_t_min,vec_t_max);
%% Debuging
if solution.problem ==1
    cons = zeros(1,N);
    for z = 1:N
        V = Dimcorrect(Data,V_result,z);
        cons(z) = violation_Not_test(Data,V);
    end
end
%Histogram_Plotting(Data,N_run,Ntest);
