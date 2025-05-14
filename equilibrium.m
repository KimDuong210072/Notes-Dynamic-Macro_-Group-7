%% File Info.

%{

    equilibrium.m
    -------------
    This code compute the equilibrium of the model.

%}

%% Solve class.

classdef equilibrium
    methods(Static)
        
        function F = obj_fun(r0,par)
            %% Find the value of r so that markets clear.

            par.r = r0; % Guess of r.

            [par,sol] = solve.firm_problem(par); % Firms.
            sol = solve.hh_problem(par,sol); % Households.
            sim = simulate.economy(par,sol);

            F = norm(sim.asup-sol.k);
            
        end
        
    end
end