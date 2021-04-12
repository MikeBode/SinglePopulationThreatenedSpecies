clear all

figure(1), clf; FS = 18;
for i = 1:4
    [D,T] = xlsread('Datasets.xlsx');
    
    %% These are the variables we're changing for the decision problem.
    % How many years has the population been observed for?
    if i == 1;      Obs = 15; 
    elseif i == 2;  Obs = 8;
    elseif i == 3;  Obs = 9;
    elseif i == 4;  Obs = 10;
    end
    ObsV(i) = Obs;
    
    %% Input the observed data
    SppName = T{1,2*i-1}
    t_raw = D(:,2*i-1); t_raw(isnan(t_raw)) = [];
    N_raw = D(:,2*i); N_raw(isnan(N_raw)) = [];
    F = find(t_raw < t_raw(1) + Obs,1,'last'); % This is the last observation before the decision point
    
    % Create an annual abundance estimate (using interpolation if necessary)
    if i < 4
        t = floor(min(t_raw)):ceil(max(t_raw));
        N = interp1(t_raw,N_raw,t,'linear');
    else
        t = t_raw;
        N = N_raw;
        F = Obs;
    end
    
    % Parameter values for the value function
    rate = log(20)/N(1); % The slope is defined by the value of the original population
    rel = 0.667; % The relative value of ex situ, in toto populations, relative to in situ populations of the same size
    w_in = [rate rel 0]; % Decay of value as abundance declines, and the relative value of in situ vs. ex situ
    w_ex = [rate rel 1]; % Decay of value as abundance declines, and the relative value of in situ vs. ex situ
    Value_function = @(x,w) w(2).^w(3).*(1 - exp(-w(1).*x));
    
    % How many years have the managers attempted the in situ action for?
    Previous_attempts = Obs-1;
    
    % Define probability that the ex situ action will be successful, based on prior experience
    p_x = 0.75;
    
    % What are the observed decline rates so far?
    ObsDecline = N(2:Obs)-N(1:Obs-1);
    MeanDecline = mean(ObsDecline);
    StdDecline = std(ObsDecline).*sqrt(Obs/(Obs-1));
    if i == 3; StdDecline = StdDecline + 0.5; end % Inflate decline rate for Christmas Island pipistrelle (no way we're that certain)
    DeclineCI = min(-1,norminv([0.025 0.975],MeanDecline,StdDecline));
    
    % Plot linear regression from the final datapoint forward
    subplot(4,1,i); hold on; box on
    t_forward = [0:30];
    y_forward_L = N(Obs) + t_forward.*DeclineCI(1);
    y_forward_H = N(Obs) + t_forward.*DeclineCI(2);
    XX = t(Obs) + [t_forward t_forward(end:-1:1)];
    YY = [y_forward_L y_forward_H(end:-1:1)];
    pp = patch(XX,YY,'r'); set(pp,'facealpha',0.2,'edgecolor','none'); PPP(1) = pp;
    plot(t,N,'k.','markersize',4)
    PPP(2) = plot(t_raw(1:F),N_raw(1:F),'ko', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot(t_raw(F+1:end),N_raw(F+1:end),'ko', 'LineWidth', 1.2, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k')
    xlim([min(t)-1 max(t)+5])
    YL = ylim; ylim([0 N(1).*1.1])
    
    anArrow = annotation('arrow') ;
    anArrow.Parent = gca;  dy = -N(1)/6;
    anArrow.Position = [t(Obs), N(Obs)-dy+N(1)/20, 0, dy] ;
    set(anArrow,'color','r','linewidth',1.5)    
    
    % Define the current population abundance when the decision is made
    N0 = N(Obs);
    Nvec = 0:N0; % This is the range of population sizes we will be forced to consider
    
    % How long will this simulation possibly run for
    Tmax = ceil(max(-N0./DeclineCI));
    if Tmax <= 0; Tmax = 25; end
    
    % Define the system states
    counter = 1;
    for Failures = 0:Tmax
        for Action = [0 1]
            for Successes = [0 1]
                
                % We won't allow ex situ action to come after success, or visa versa
                StateVec(counter,:) = [Failures Action Successes];
                
                % What's the expected abundance in this state?
                ElapsedTime(counter) = sum(StateVec(counter,:));
                
                % What's the pdf of the abundance at this point in time?
                UCI(counter,1) = max(0,N0 + DeclineCI(2).*ElapsedTime(counter));
                LCI(counter,1) = max(0,N0 + DeclineCI(1).*ElapsedTime(counter));
                mu(counter,1)  = max(0,N0 + MeanDecline.*ElapsedTime(counter));
                sigma(counter,1) = max(1e-1,abs(UCI(counter,1)-LCI(counter,1))/2);
                Abundance(counter,:) = normpdf(Nvec,mu(counter,1),sigma(counter,1));
                Abundance(counter,:) = Abundance(counter,:)./sum(Abundance(counter,:));
                ExpectedN(counter,1) = sum(Abundance(counter,:).*Nvec);
                
                % Iterate the counter
                counter = counter + 1;
            end
        end
    end
    NumStates = size(StateVec,1);
    
    % Define the value vector
    Value = zeros(NumStates,Tmax);
    for s = 1:NumStates
        
        if StateVec(s,3) == 1 & StateVec(s,2) == 1  % Transition implies that we succeeded in situ and ex situ, which is impossible
            Value(s,end) = 0;
            continue
        end
        
        if StateVec(s,3) == 1 % We succeeded in situ
            Value(s,end) = sum(Abundance(s,:).*Value_function(Nvec,w_in));
            continue
        end
        
        if StateVec(s,2) == 1 % We succeeded ex situ
            Value(s,end) = sum(Abundance(s,:).*Value_function(Nvec,w_ex));
            continue
        end
        
        % We're still waiting for something to work
        Value(s,end) = 0;
    end
    % Define the final state as extinction
    Value(end,end) = 0;
    ExpectedN(end) = -inf;
    
    %% =============== In situ actions ===============
    T_in = zeros(NumStates); %
    % StateVec(counter,:) = [In-situ-Failures     Ex-situ-Action     In-situ-Successes];
    for sI = 1:NumStates
        for sF = 1:NumStates
            
            % If we took an action and it worked
            if (StateVec(sI,3) == 0 & StateVec(sF,3) == 1) ... % Success went from zero to one
                    & (StateVec(sF,1) == StateVec(sI,1)) ... % Failures didn't increase
                    & (StateVec(sF,2) == 0 & StateVec(sI,2) == 0) % Ex situ was never taken.
                T_in(sI,sF) = 1/(StateVec(sI,1)+Previous_attempts + 2);
            end
            
            % If we took an action and it failed
            if (StateVec(sI,3) == 0 & StateVec(sF,3) == 0) ... % Success stayed at zero
                    & (StateVec(sF,1) == StateVec(sI,1)+1) ... % Failures increased by one
                    & (StateVec(sF,2) == 0 & StateVec(sI,2) == 0) % Ex situ was never taken.
                T_in(sI,sF) = 1 - 1/(StateVec(sI,1)+Previous_attempts + 2);
            end
            
            if (StateVec(sI,2) == 1 | StateVec(sI,3) == 1) & (sI == sF) % Actions have already been taken
                T_in(sI,sF) = 1;
            end
            
        end
    end
    
    % Deal with one final state (there's nowhere for {\tau,0,0} to go if the action fails)
    T_in(NumStates-3,NumStates) = 1 - T_in(NumStates-3,NumStates-2);
    
    %% =============== Ex situ, in toto actions ===============
    T_ex = zeros(NumStates);
    % StateVec(counter,:) = [In-situ-Failures     Ex-situ-Action      In-situ-Successes];
    for sI = 1:NumStates
        for sF = 1:NumStates
            
            % If we took ESIT action and it worked
            if (StateVec(sI,2) == 0 & StateVec(sF,2) == 1) ... % Ex situ hasn't already been done, and success went from zero to one
                    & (StateVec(sF,1) == StateVec(sI,1)) ... % Failures didn't increase
                    & (StateVec(sF,3) == 0 & StateVec(sI,3) == 0) % In situ hasn't previously worked
                T_ex(sI,sF) = p_x;
            end
            
            % If we took ESIT action and it failed
            if (StateVec(sI,2) == 0 & StateVec(sI,3) == 0) ... % Ex situ and in situ haven't been done
                    & sF == 4*(Tmax+1) % And we're going straight to the end
                T_ex(sI,sF) = 1-p_x;
            end
            
            if (StateVec(sI,2) == 1 | StateVec(sI,3) == 1) & (sI == sF) % Actions have already been taken
                T_ex(sI,sF) = 1;
            end
            
        end
    end
    
    % Check that these are approximately probability matrices
    if min(sum(T_in,2)) < 0.999 | min(sum(T_ex,2)) < 0.999
        disp('Error in transition matrix')
    end
    
    %% Apply SDP
    for nbs = 0:Tmax-2
        Options(:,1) = T_in*Value(:,Tmax-nbs); % Continued in situ action
        Options(:,2) = T_ex*Value(:,Tmax-nbs); % Captive breeding
        [Value(:,Tmax-nbs-1),OptimalAction(:,Tmax-nbs-1)] = max(Options,[],2);
    end
    % Relevant states have neither success or failure yet
    RelSt = find(sum(StateVec(:,2:3),2)==0);
    
    if i == 1
        xlim([1990 2022])
    elseif i == 2
        xlim([1982 1998])
    elseif i == 3
        xlim([1993 2010])
    elseif i == 4
        xlim([1978 1996])
    end
    
    if i == 1
        L = legend(PPP(2:-1:1),'Observations','Prediction 95\% CIs');
        set(L,'location','northeast','fontsize',FS-2,'interpreter','latex','box','off')
    end
    if i == 4
        xlabel('Survey year','interpreter','latex','fontsize',FS)
    end
    ylabel('Abundance','interpreter','latex','fontsize',FS)
    
    XL = xlim; YL = ylim; G = [0.04 0.1];
    text(XL(1)+G(1)*(XL(2)-XL(1)),YL(1)+G(2)*(YL(2)-YL(1)),SppName,'fontweight','bold','interpreter','latex','fontsize',FS-3,'color',[0 0 0.5])    
    
    
    clearvars -except i ha FS ObsV
end








