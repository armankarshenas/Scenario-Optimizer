function [] = Plot_Optimal(Data,t,vec_t_min,vec_t_max)
%PLOT_OPTIMAL plots the optimal solution on top of all the test scenarios'
%responses
X_result = Data.Xr;
x_o = Data.xi;
% Zone 1
subplot(3,1,1);
scatter(t,[x_o(1,1) X_result(1,:)],'r'); title('Temperature of Zone 1'); xlabel('Time / Hours'); ylabel('Temperature / C');hold on;
plot(t,[x_o(1,1) X_result(1,:)],'r'); hold on;
if Data.Summer ==1
    plot(t,vec_t_max,'r--'); hold on
elseif Data.Summer ==0
    plot(t,vec_t_min,'r--'); hold on
else
    fprintf('Season must be selected');
end
% Zone 2
subplot(3,1,2);
scatter(t,[x_o(2,1) X_result(2,:)],'r'); title('Temperature of Zone 2'); xlabel('Time / Hours'); ylabel('Temperature / C');hold on;
plot(t,[x_o(2,1) X_result(2,:)],'r'); hold on;
if Data.Summer ==1
    plot(t,vec_t_max,'r--'); hold on
elseif Data.Summer ==0
    plot(t,vec_t_min,'r--'); hold on
else
    fprintf('Season must be selected');
end

% Zone 3
subplot(3,1,3);
scatter(t,[x_o(3,1) X_result(3,:)],'r'); title('Temperature of Zone 3'); xlabel('Time / Hours'); ylabel('Temperature / C');hold on;
plot(t,[x_o(3,1) X_result(3,:)],'r'); hold on;
if Data.Summer ==1
    plot(t,vec_t_max,'r--'); hold on
elseif Data.Summer ==0
    plot(t,vec_t_min,'r--'); hold on
else
    fprintf('Season must be selected');
end


end

