%% File Info.

%{
    solve.m
    -------
    This code solves the model.
%}

%% Solve class.

classdef solvewithtax
    methods(Static)

        function sol = hh_problem(par,sol)            
            beta = par.beta;
            agrid = par.agrid;
            alen = par.alen;
            zgrid = par.zgrid;
            zlen = par.zlen;
            pmat = par.pmat;
            r = par.r;
            w = par.w;
            tau = par.tau;
            phi = 0; % No borrowing

            v1 = nan(alen,zlen,2);
            a1 = nan(alen,zlen);
            c1 = nan(alen,zlen);
            c2 = nan(alen,zlen); % Young consumption
            c3 = nan(alen,zlen); % Old consumption


            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;

            fprintf('------------Beginning Value Function Iteration.------------\n\n')

            c0 = (1 + r) * agrid + w .* zgrid;
            v0 = modelwithtax.utility(c0, par) ./ (1 - beta);
            
            while diff > crit && iter < maxiter %Cần sửa: change to for loop, vd như for age= 1-10 or smthing (j = 10)
                for i = 1:alen
                    if agrid(i) >= phi
                        for j = 1:zlen
                            income = (1 - tau) * w * par.e(2) * zgrid(j); % middle-aged only
                            T = tau * w * mean(zgrid); % Transfer to young
                            asset_choices = agrid;
                            c = (1 + r) * agrid(i) + income + T - asset_choices;
        
                            ev = v0 * pmat(j,:)';
                            vall = modelwithtax.utility(c, par) + beta * ev;
                            vall(c <= 0) = -inf;
        
                            [vmax, ind] = max(vall);
        
                            v1(i,j) = vmax;
                            c1(i,j) = c(ind);
                            a1(i,j) = asset_choices(ind);
        
                            c2(i,j) = T; % Young consumption
                            c3(i,j) = a1(i,j); % Old consumption
                        end
                    end
                end
        
                diff = norm(v1 - v0);
                v0 = v1;
                iter = iter + 1;
        
                if mod(iter, 25) == 0
                    fprintf('Iteration: %d.\n', iter)
                end
            end

            fprintf('\nConverged in %d iterations.\n\n', iter)
            fprintf('------------End of Value Function Iteration.------------\n')

            sol.a = a1;
            sol.c = c1;
            sol.c2 = c2;     % Young consumption
            sol.c3 = c3;     % Old consumption
            sol.v = v1;
        end

        function [par,sol] = firm_problem(par)
            delta = par.delta;
            alpha = par.alpha;
            r = par.r;
            k = ((r + delta) / alpha) ^ (1 / (alpha - 1));

            sol = struct();
            sol.k = k;
            par.w = (1 - alpha) * k ^ alpha;
        end
    end
end

%UBI: In the budget constraint, 
%Consider linear taxes
%Welfare effect, compute average utility of consumption. 
%