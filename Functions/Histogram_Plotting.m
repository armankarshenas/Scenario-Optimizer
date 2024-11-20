function [viol_vec] = Histogram_Plotting(Data,N_run,Ntest)
%HISTOGRAM_PLOTTING Plots a histogram for the distribution in scenario
%violation
viol_vec = zeros(1,N_run);
for i = 1 : N_run
    viol_vec(i) = violation(Data,Ntest,0);
end
%figure; title('Violation distribution');
%histogram(viol_vec/Ntest);

end

