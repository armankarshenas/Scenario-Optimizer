%% Scenario approach in MPC for energy management in a building
%% Initation
close all;
clear;
clc;
addpath(genpath('/Users/Karshenas/Arman/Research + Paper/2019-2020/EUROP 2019/tbxmanager/toolboxes'));

% Sampling time in hours
Sampling = 0.25;
Data.Para.Sampling = Sampling;

% Starting the clock for execution time
tic;

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

% The input variables Gamma and Theta for the feedforward and the feedback
% gain respectively
Gamma = sdpvar(M*u_length,1);
Theta = sdpvar(M*u_length,M*v_length);
Data.Gamma = Gamma;
Data.Theta = Theta;

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

% Matrices to store scenarios
X_result = zeros(x_length,M);
V_result = zeros(v_length,M,N);
U_result = zeros(u_length,M);

% Testing Scenarios parametes
Ntest = 3000;
Nrun = 100;

% The uncertainty vector parameters
lambda  = 3;
Amb_temp = Data.Para.Amb;
Global_radiation = 200;
Data.Para.Lambda = lambda;
Data.Para.GR = Global_radiation;

% index vectors used in complexity computation
indexcountermax= zeros(2,M);
indexcountermin= zeros(2,M);
%% Incremental Parameters
% Incremental Scenario Approach Parameters
complexity = 1;

% Vector N
n_s = Numberscenario(Data);
n_s_oo = 0;

% The expected value and the second moment of uncertainty used in the cost
expected =Data.Para.Lambda*ones(M*v_length,1)*9.324;
secondmoment = (Data.Para.Lambda*9.324+(Data.Para.Lambda*9.324)^2)*ones(M*v_length,M*v_length);

% Cost function with feedback
Cost = Gamma'*R*Gamma + 2*Gamma'*R*Theta*expected+trace(Theta'*R*Theta*secondmoment);
%% Constraints on Theta
Constraints = [];
zero_theta = zeros(u_length,v_length);
% Setting Theta to a lower triangular matrix
for i=0:M-1
    for j=i:M-1
        Constraints = [ Constraints Theta(i*u_length+1:i*u_length + u_length,j*v_length+1:j*v_length+v_length) == zero_theta];
    end
end
%% Incremental Algorithm and Solver
j=0;
while (complexity > j) && (j< M*u_length)
    % Getting Scenarios if j ==0
    if j==0
        maxv=zeros(M,1);
        minv=60*ones(M,1);
        
        for k=1:n_s(j+1)
            [V,v_local] = ScenarioGene(Data);
            V_result(:,:,k) =v_local;
        end
        
        for i=1:M
            for k=1:n_s(j+1)
                V = Dimcorrect(Data,V_result,k);
                if V((i-1)*v_length+1,1) > maxv(i)
                    maxv(i) = V((i-1)*v_length+1,1);
                    indexcountermax(1,i) =k;
                end
                if V((i-1)*u_length+1,1) < minv(i)
                    minv(i) = V((i-1)*u_length+1,1);
                    indexcountermin(1,i) =k;
                end
                
            end
        end
        
        Vmax = V;
        for i=0:M-1
            Vmax(i*v_length+1,1) = maxv(i+1);
        end
        
        Vmin = V;
        for i=0:M-1
            Vmin(i*u_length+1,1) = minv(i+1);
        end
        if Summer ==1
            % The constraints introduced by scenarios
            Constraints = [Constraints, F*x_o+G*Gamma+(H+G*Theta)*Vmax <= g_t_max ];
            
        elseif Summer ==0
            % The constraints introduced by scenarios
            Constraints =[Constraints, F*x_o+G*Gamma+(H+G*Theta)*Vmin >= g_t_min];
        else
            fprintf('Season must be selected');
        end
    else
        % The uncertainty vector gets new scenarios now that j != 0
        maxv=zeros(M,1);
        for k=n_s(j+1):n_s(j+2)
            [V,v_local] = ScenarioGene(Data);
            V_result(:,:,k) =v_local;
        end
        minv=V(M*v_length-2,1)*ones(M,1);
        for i=1:M
            for k=n_s(j+1):n_s(j+2)
                V = Dimcorrect(Data,V_result,k);
                if V((i-1)*u_length+1,1) > maxv(i)
                    maxv(i) = V((i-1)*u_length+1,1);
                    indexcountermax(1,i) =k;
                end
                if V((i-1)*u_length+1,1) < minv(i)
                    minv(i) = V((i-1)*u_length+1,1);
                    indexcountermin(1,i) =k;
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
            % The constraints introduced by scenarios
            Constraints = [Constraints, F*x_o+G*Gamma+(H+G*Theta)*Vmax <= g_t_max ];
            
        elseif Summer ==0
            % The constraints introduced by scenarios
            Constraints =[Constraints, F*x_o+G*Gamma+(H+G*Theta)*Vmin >= g_t_min];
        else
            fprintf('Season must be selected');
        end
    end
    % Input Constraints
    if Summer ==1
        [Constraintsinput,Data] = Input_Constraints_feedback(Data,Vmax);
    elseif Summer ==0
        [Constraintsinput,Data] = Input_Constraints_feedback(Data,Vmin);
    end
    
    Constraints = [Constraints , Constraintsinput];
    
    % Optimisation
    options = sdpsettings('solver','quadprog','debug',1);
    solution = optimize(Constraints,Cost,options);
    U = value(Gamma)+value(Theta)*Vmax;
    Data.U= U;
    Activescenario = zeros(1,M*x_length);
    
    % Finding the Complexity of the solution
    complexityini =0;
    X_max = F*x_o+G*U+H*Vmax;
    X_min = F*x_o+G*U+H*Vmin;
    UB = Data.UB;
    LB = Data.LB;
    for p=0: M*x_length-1
        % Checking for active constraints within a threshhold of the bounds
        
        if ((X_max(p+1,1) < g_t_max(p+1,1)+1e-3) && (X_max(p+1,1) > g_t_max(p+1,1)-1e-3)) || ((U(floor(u_length/x_length)*p+1,1) >=LB(floor(u_length/x_length)*p+1,1)-1e-3) && (U(floor(u_length/x_length)*p+1,1) <=LB(floor(u_length/x_length)*p+1,1)+1e-3)) || ((U(floor(u_length/x_length)*p+1,1) >=UB(floor(u_length/x_length)*p+1,1)-1e-3) && (U(floor(u_length/x_length)*p+1,1) <=UB(floor(u_length/x_length)*p+1,1)+1e-3))
            complexityini = complexityini +1;
            Activescenario(p+1) = indexcountermax(1,floor(p/x_length)+1);
        end
        
    end
    % making sure that the same scenario is not counted more than once
    
    complexity = length(unique(Activescenario))-1;
    
    if (complexity <= j)
        U = value(Gamma);
        break;
    else
        j = j +1;
    end
    
end



%% Solution treatment
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
%% Plotting States, Inputs and Uncertainties
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
viol_vec = Histogram_Plotting(Data,Nrun,Ntest);
meanvec = mean(viol_vec);
Varvec = var(viol_vec);
f = @(x) (1/sqrt(2*pi*Varvec))*(exp(-0.5*((x-meanvec).^2/Varvec)));
y = f(viol_vec);
histogram(viol_vec); hold on;
scatter(viol_vec,y,'r');
toc;