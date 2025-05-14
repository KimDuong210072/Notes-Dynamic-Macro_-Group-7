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
            phi = 0; % No borrowing

            v1 = nan(alen,zlen);
            a1 = nan(alen,zlen);
            c1 = nan(alen,zlen);

            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;

            fprintf('------------Beginning Value Function Iteration.------------\n\n')

            c0 = (1 + r) * agrid + w .* zgrid;
            v0 = model.utility(c0, par) ./ (1 - beta);
            
            while diff > crit && iter < maxiter
                for p = 1:alen
                    if agrid(p) >= phi
                        for j = 1:zlen
                            income_after_tax = (1 - par.tau) * w * zgrid(j);
                            T = par.tau * w * mean(zgrid); % Lump-sum transfer
                            c = (1 + r) * agrid(p) + income_after_tax + T - agrid;

                            ev = v0 * pmat(j,:)';
                            vall = model.utility(c, par) + beta * ev;
                            vall(c <= 0) = -inf;
                            [vmax, ind] = max(vall);

                            v1(p,j) = vmax;
                            c1(p,j) = c(ind);
                            a1(p,j) = agrid(ind);
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
