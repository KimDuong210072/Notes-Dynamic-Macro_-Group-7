%% File Info.

%{

    equilibrium.m
    -------------
    This code compute the equilibrium of the model.

%}

%% Solve class.

classdef lequilibrium
    methods(Static)
        
        function F = obj_fun(r0,par)
            %% Find the value of r so that markets clear.

            par.r = r0; % Guess of r.

            [par,sol] = lsolvewithtax.firm_problem(par); % Firms.
            sol = lsolvewithtax.hh_problem(par,sol); % Households.
            sim = lsimulate.economy(par,sol);

            F = norm(sim.asup-sol.k);
            
        end
        
    end
end