function [F,G,H,Q,R,Z,Data] = MatrixGene(Data)
%INITATION This function Generates Matrices required for MPC

M = Data.Para.Horizon;
%% Dynamic Matrices
[A,Bu,Bv,Bvu,Bxu , Building] = dynamic(Data.Para.Sampling);
x_length = length(Building.building_model.identifiers.x);
v_length = length(Building.building_model.identifiers.v);
u_length = length(Building.building_model.identifiers.u);
Data.Building = Building;
%% Define F,G,H  x_{+} = Fx + Gu +Hv
% Matrix F which is Mn*n
F = zeros(M*x_length,x_length);
for i=0:M-1
    F(i*x_length+1:i*x_length+x_length,:) = A^(i+1);
end
Data.F = F;
% Matrix G which is Mn * Mm
G = zeros(M*x_length,M*u_length);
for i =0 : (M-1)
    G(i*x_length+1:i*x_length+x_length,i*u_length+1 : i*u_length+u_length) = Bu;
    for j=0 :i-1
        G(i*x_length+1:i*x_length+x_length,j*u_length+1 : j*u_length+u_length) = A^(i-1)*Bu;
    end
end
Data.G = G;
% Matrix H which is Mn * Ml
H = zeros(M*x_length, M*v_length);
for i =0 : (M-1)
    H(i*x_length+1:i*x_length+x_length,i*v_length+1 : i*v_length+v_length) = Bv;
    for j=0 :i-1
        H(i*x_length+1:i*x_length+x_length,j*v_length+1 : j*v_length+v_length) = A^(i-1)*Bv;
    end
end
Data.H = H;
% Matrix Q which is Mn * Mn
Q_p = zeros(x_length,x_length,M);
Q = zeros(M*x_length,M*x_length);
for i=0 : M-1
    Q(i*x_length+1 :i*x_length+x_length,i*x_length+1:i*x_length+x_length) = Q_p(:,:,i+1);
end
Data.Q = Q;
% Matrix R which is Mm * Mm
R_p = eye(u_length,u_length);
R = zeros(M*u_length,M*u_length);
for i=0 :(M-1)
    R(i*u_length+1:i*u_length+u_length,i*u_length+1:i*u_length+u_length) = R_p;
end
Data.R = R;
% Matrix Z for state constraints
Z = eye(M*x_length,M*x_length);

end

