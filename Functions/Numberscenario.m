function n_s = Numberscenario(Data)
% This function generates a vector for the number of scenarios required at
% every ietration of the incremental scenario approach

% The maximum complexity
d = Data.Para.Horizon*Data.Para.ul;
M_j = zeros(1,d+1);
n_s = zeros(1,d+1);
% Calculating Mjs
for j=0:d
    M_j(j+1) = (2/Data.Para.Eps)*(j-1+log(1/Data.Para.Beta));
    sum = 0;
    
    for m=j:floor(M_j(j+1))
        sum = sum + nchoosek(m,j)*(Data.Para.Beta*((1-Data.Para.Eps)^(m-j)))/(((d+1)*(M_j(j+1)+1)));
        clc;
    end
    % Calculating Njs
    n_s(j+1) = ceil((2/Data.Para.Eps)*log(1/sum)+2*j+(2*j/Data.Para.Eps)*log(2/Data.Para.Eps));
    
end
end