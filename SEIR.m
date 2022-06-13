classdef SEIR 
    % SEIR Mutli-group epidemic model.
    
    properties
        n       % number of groups
        beta    % transmission rate 
        gamma   % 1 / avg. latent infection time
        delta   % 1 / avg. infectious time
        s0      % Equilibrium population
        A       % Contact rate matrix
    end
    
    methods
        
        function model = SEIR(beta, gamma, delta, s0, A)
            % SEIR Create a new SEIR model.
            % Inputs:
            %   beta    : transmission rates in each group
            %   gamma   : E -> I transition rates
            %   delta   : removal rate of infected individuals
            %   s0      : equilibrium populations
            %   A       : inter-group contact rate matrix
            
            model.n = size(A, 1);
            model.beta = beta; 
            model.gamma = gamma;
            model.delta = delta;
            model.s0 = s0;
            model.A = A; 
        end
               
        function r0 = R0(model)
            % R0 Compute the model's basic reproduction number.
            
            [F, V] = model.NGM();
            r0 = max(abs(eig(F / V)));
        end 
        
        function alpha = abscissa(model)
            % ABSCISSA Compute the spectral abscissa of the model.
            
            [F, V] = model.NGM();
            alpha = max(real(eig(F + V))); 
        end
        
        function [F, V] = NGM(model)
            % NGM Compute the next generation matrix components of the
            % model.
            % Outputs: 
            %   F   : Matrix of new infection rates 
            %   V   : Matrix of transition rates not corresponding to new
            %         infections 
            
            F = [zeros(model.n) diag(model.beta .* model.s0) * model.A;
                zeros(model.n, 2*model.n)];
            V = [-diag(model.gamma) zeros(model.n); 
                diag(model.gamma) -diag(model.delta)];   
        end
        
        function infx = totalInfections(model, x0)
            % TOTALINFECTIONS Compute the total cumulative number of 
            % infections in the t -> \infty limit (not counting initial 
            % infections).
            % Inputs:
            %   x0  : initial condition (s0, e0, z0, r0)
            % Outputs:
            %   infx: number of infections
            
            e0 = x0(1 + model.n : 2 * model.n);
            z0 = x0(1 + 2 * model.n : 3 * model.n);
            [F, V] = model.NGM();
            q = (F + V) \ [e0; z0];
            infx = -F(1 : model.n, :) * q; 
        end
        
        function [t, s, e, z, r] = simulate(model, x0, tf) 
            % SIMULATE Simulate the nonlinear dynamics.
            % Inputs:
            %   x0  : initial condition (s0, e0, z0, r0)
            %   tf  : time duration of simulation
            % Outputs:
            %   t   : simulation times
            %   s   : simulated numbers of susceptible individuals
            %   e   : simulated numbers of exposed individuals
            %   z   : simulated numbers of infected individuals
            %   r   : simulated numbers of removed individuals
            
            fn = @(t, x) model.fn(x);
            [t, x] = ode45(fn, [0 tf], x0); 
            s = x(:, 1 : model.n);
            e = x(:, 1 + model.n : 2 * model.n);
            z = x(:, 1 + 2 * model.n : 3 * model.n);
            r = x(:, 1 + 3 * model.n : 4 * model.n); 
        end      
        
        function dx = fn(model, x)
             
            s = x(1:model.n);
            e = x(1+model.n:2*model.n);
            z = x(1+2*model.n:3*model.n); 
            
            newIfx = diag(model.beta) * diag(s) * (model.A * z);
            ds = -newIfx;
            de = newIfx - model.gamma .* e;
            dz = model.gamma .* e - model.delta .* z;
            dr = model.delta .* z;
            dx = [ds; de; dz; dr];
        end
        
    end
end

