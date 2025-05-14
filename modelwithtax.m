%% File Info.

%{
    model.m
    -------
    This code sets up the model.
%}

%% Model class.

classdef modelwithtax
    methods(Static)

        function par = setup()            
            par = struct();

            %% Preferences.
            par.beta = 0.94;
            par.sigma = 2.00;
            assert(par.beta > 0 && par.beta < 1.00)
            assert(par.sigma > 0)

            %% Technology.
            par.alpha = 0.33;
            par.delta = 0.03;
            assert(par.alpha > 0 && par.alpha < 1.00)
            assert(par.delta >= 0 && par.delta <= 1.00)

            par.sigma_eps = 0.07;
            par.rho = 0.85;
            par.mu = 0.0;
            assert(par.sigma_eps > 0)
            assert(abs(par.rho) < 1)

            %% Prices.
            par.r = 0.06;
            par.omega = 0.2;
            par.tau = 0.2; % Income tax rate: 20%

            %% Simulation parameters.
            par.seed = 2025;
            par.T = 10000;
            par.N = 3000;
        end
        
        function par = gen_grids(par)
            par.alen = 500;
            par.amax = 300;
            par.amin = 0;
            assert(par.alen > 5)
            assert(par.amax > par.amin)
            par.agrid = linspace(par.amin,par.amax,par.alen)';

            par.zlen = 7;
            par.m = 3;
            assert(par.zlen > 3)
            assert(par.m > 0)

            [zgrid,pmat] = model.tauchen(par.mu,par.rho,par.sigma_eps,par.zlen,par.m);
            par.zgrid = exp(zgrid);
            par.pmat = pmat;
        end

        function [y,pi] = tauchen(mu,rho,sigma,N,m)
            ar_mean = mu/(1-rho);
            ar_sd = sigma/((1-rho^2)^(1/2));
            y1 = ar_mean - (m * ar_sd);
            yn = ar_mean + (m * ar_sd);
            y = linspace(y1, yn, N);
            d = y(2) - y(1);
            ymatk = repmat(y, N, 1);
            ymatj = mu + rho * ymatk';
            pi = normcdf(ymatk, ymatj - (d/2), sigma) - normcdf(ymatk, ymatj + (d/2), sigma);
            pi(:,1) = normcdf(y(1), mu + rho * y - (d/2), sigma);
            pi(:,N) = 1 - normcdf(y(N), mu + rho * y + (d/2), sigma);
        end

        function u = utility(c,par)
            if par.sigma == 1
                u = log(c);
            else
                u = (c.^(1-par.sigma))./(1-par.sigma);
            end
        end
    end
end
