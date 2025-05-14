%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions.

%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_dist_pe(sim)
            %% Plot density function (partial equilibrium).
            
            figure(1)
            
            histogram(sim.asim,20)
            xlabel({'$a$'},'Interpreter','latex')
            ylabel({'Frequency'},'Interpreter','latex') 
            title('Distribution of Wealth')

        end
            
        function [] = plot_dist_ge(sim_ge)
            %% Plot density function (general equilibrium).
            
            figure(2)
            
            histogram(sim_ge.asim,20)
            xlabel({'$a$'},'Interpreter','latex')
            ylabel({'Frequency'},'Interpreter','latex') 
            title('Distribution of Wealth')

        end
            
        function [] = cfun(par,sol)
            
            %% Plot consumption policy function.
            
            figure(3)
            
            plot(par.agrid,sol.c)
            xlabel({'$a_{t}$'},'Interpreter','latex')
            ylabel({'$c_{t}$'},'Interpreter','latex') 
            title('Consumption Policy Function')
            
        end
        
    end
end