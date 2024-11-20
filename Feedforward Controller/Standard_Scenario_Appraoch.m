%% Scenario approach in MPC for energy management in a building  
%% Initation and Dynamic constraint
close all;
clear;
clc;
addpath(genpath('/Users/Karshenas/Arman/Research + Paper/2019-2020/EUROP 2019/tbxmanager/toolboxes'));
% Sampling time in hours
Sampling = 0.25;
Data.Para.Sampling = Sampling;

% MPC Horizon and scenarios
Horizon_time = 12;
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

% The input variable Gamma
Gamma = sdpvar(M*u_length,1);
Data.Gamma = Gamma;

% Summer/Winter:
Summer =1;
Data = Season(Data,Summer);
Data.Summer = Summer;

% The Initial Condition
x_o = Data.xi;

% Optimal Temperature Range
t_min = Data.Para.Tmin;
t_max = Data.Para.Tmax;

% The temperature range vectors used for comfort restrictions
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

%% Input Constraints
% Constraints on the input
Constraints =[];
Constraints = Input_Constraints(Data);

%% Scenarios and Constrainst 
% finding the worst possible uncertainty vectors from the scenario sample 
maxv=zeros(M,1);

for k=1:N
    [V,v_local] = ScenarioGene(Data);
    V_result(:,:,k) =v_local;
end
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
end
if Summer ==1 
    
Constraints = [Constraints, F*x_o+G*Gamma+H*Vmax <= g_t_max ];
elseif Summer ==0
    Constraints =[Constraints, F*x_o+G*Gamma+H*Vmin >= g_t_min];
else 
    fprintf('Season must be selected');
end

% Cost Function
Cost = Gamma'*R*Gamma;

%Optimisation
options = sdpsettings('solver','quadprog','debug',1);
solution = optimize(Constraints,Cost,options);
U = value(Gamma);
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
%% Histogram Plotting 
Histogram_Plotting(Data,N_run,Ntest);
