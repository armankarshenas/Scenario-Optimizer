function violcount = violation(Data,Ntest,draw)
% This function tests the solution obtaiend with a new set of scenarios and
% counts the number of occasaions that the optimal solution violates the
% constraints
violcount = 0;
X_result = zeros(Data.Para.xl,Data.Para.Horizon);
for k = 1 :Ntest
    % getting a new realisation
    [V,v_local] = ScenarioGene(Data);
    X = Data.F*Data.xi+Data.G*Data.U+Data.H*V;
    test1 = Data.Tmax;
    test2 = Data.Tmin;
    viol =0;
    for i=0 : Data.Para.Horizon-1
        if Data.Summer ==1
            if (X(i*Data.Para.xl+1,1) > test1(1,1)) || (X(i*Data.Para.xl+2,1) > test1(2,1)) || (X(i*Data.Para.xl+3,1) > test1(3,1))
                if (viol ==0)
                    violcount = violcount +1;
                    viol=1;
                end
            end
        end
        if Data.Summer==0
            if (X(i*Data.Para.xl+1,1) < test2(1,1)) || (X(i*Data.Para.xl+2,1) < test2(2,1)) || (X(i*Data.Para.xl+3,1) < test2(3,1))
                if (viol ==0)
                    violcount = violcount +1;
                    viol=1;
                end
            end
        end
        for j=0:Data.Para.Horizon-1
            X_result(:,j+1) = X(j*Data.Para.xl+1:j*Data.Para.xl+Data.Para.xl,1);
        end
    end
    
    % if draw ==1, the function plots the response of the system to each of
    % the realisations of the uncertainty
    if draw == 1
        t = 0 : Data.Para.Sampling : Data.Para.Horizon*Data.Para.Sampling;
        subplot(3,1,1);hold on;
        plot(t,[Data.xi(1,1) X_result(1,:)],'color',[209/255 200/255 200/255],'LineWidth',3); hold on;
        
        subplot(3,1,2);hold on
        plot(t,[Data.xi(2,1) X_result(2,:)],'color',[209/255 200/255 200/255],'LineWidth',3); hold on;
        
        subplot(3,1,3);hold on;
        plot(t,[Data.xi(3,1) X_result(3,:)],'color',[209/255 200/255 200/255],'LineWidth',3); hold on;
        
    end
end
end
