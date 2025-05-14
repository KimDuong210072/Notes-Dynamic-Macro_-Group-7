%% File Info.

%{

    simulate.m
    ----------
    This code simulates the model.

%}

%% Simulate class.

classdef simulate
    methods(Static)
        %% Simulate the model. 
        
        function sim = economy(par,sol)            
            %% Set up.
            
            agrid = par.agrid; % Assets today (state variable).
            zgrid = par.zgrid; % Productivity.
            pmat = par.pmat; % Transition matrix.

            apol = sol.a; % Policy function for capital.
            cpol = sol.c; % Policy function for consumption.

            T = par.T+10; % Time periods.
            N = par.N; % People.

            zsim = nan(T,N); % Container for simulated income.
            asim = nan(T,N); % Container for simulated savings.
            csim = nan(T,N); % Container for simulated consumption.
            usim = nan(T,N); % Container for simulated utility.
            
            %% Begin simulation.
            
            rng(par.seed);

            pmat0 = pmat^100; % Stationary distirbution.
            cmat = cumsum(pmat,2); % CDF matrix.

            z0_ind = randsample(par.zlen,N,true,pmat0(1,:))'; % Index for initial income.
            a0_ind = randsample(par.alen,N,true)'; % Index for initial wealth.

            for i = 1:N % Person loop.
                zsim(1,i) = zgrid(z0_ind(i)); % Productivity in period 0.
                csim(1,i) = cpol(a0_ind(i),z0_ind(i)); % Consumption in period 0 given a0.
                asim(1,i) = apol(a0_ind(i),z0_ind(i)); % Savings for period 1 given a0.
                usim(1,i) = model.utility(csim(1,i),par); % Utility in period 0 given a0.
            end

            %% Simulate endogenous variables.

            for j = 2:T % Time loop.
                for i = 1:N % Person loop.
                    at_ind = find(asim(j-1,i)==agrid); % Savings choice in the previous period is the state today. Find where the latter is on the grid.
                    zsim(j,i) = zgrid(z0_ind(i)); % Productivity in period 0.
                    csim(j,i) = cpol(at_ind,z0_ind(i)); % Consumption in period t.
                    asim(j,i) = apol(at_ind,z0_ind(i)); % Savings for period t+1.
                    usim(j,i) = model.utility(csim(j,i),par); % Utility in period t.
                    z1_ind = find(rand<=cmat(z0_ind(i),:)); % Draw income shock for next period.
                    z0_ind(i) = z1_ind(1);
                end

                if j >= 11 && norm(mean(asim(j-9:j,:),'all')-mean(asim(j-10:j-1,:),'all')) < 0.0001
                    break
                end

            end

            sim = struct();
            
            sim.zsim = zsim(j-9:j,:); % Simulated productivity.
            sim.asim = asim(j-9:j,:); % Simulated savings.
            sim.csim = csim(j-9:j,:); % Simulated consumption.
            sim.usim = usim(j-9:j,:); % Simulated utility.
            sim.asup = mean(asim(j-9:j,:),'all'); % Simulated savings.
             
        end
        
    end
end