classdef SEIRAllocator
    % SEIRALLOCATOR Allocates resources (vaccines, antidotes, NPIs) in a 
    % given pre-intervention SEIR model. There are 2 main methods:
    %   allocatePharmaR0        : budget-constrained allocation of vaccines
    %                             and antidotes to minimize R0
    %   allocatePharmaAbscissa  : budget-constrained allocation of vaccines
    %                             and antidotes to minimize the spectral
    %                             absicssa
    
    properties
        
        % Allocator parameters
        model       % Pre-intervention SEIR model
        betaRange   % [betaMin, betaMax] for each node
        deltaRange  % [deltaMin, deltaMax] for each node 
        delOffset
        
        % Pre-computed cost function parameters
        fSlope      % Coefs of 1./beta in vax cost fn
        gSlope      % Coefs of 1./(1-delta) in antidote cost fn 
        fOffset     % Vax costs when beta are maximal
        gOffset     % Antidote costs when delta are minimal 
    end
    
    methods
        
        function alloc = SEIRAllocator(model, betRange, delRange, delOffset) 
            % SEIRALLOCATOR Create a new resource allocator for the SEIR
            % model.
            % Inputs:
            %   model       : pre-intervention SEIR model
            %   betRange    : lower bounds betaRange(:,1) and upper bounds
            %                 betaRange(:,2) on beta
            %   delRange    : lower bounds deltaRange(;,1) and upper bounds
            %                 deltaRange(:,2) on delta 
            %   delOffset   : upper bound on the values of delta
            
            % Check inputs 
            if delOffset <= max(delRange(:,2))
                error([...
                    'Delta offset ' num2str(delOffset) ...
                    ' is too small to transcribe the optimization' ...
                    ' as a geometric program. Try a value of' ...
                    ' delOffset that is greater than ' ...
                    num2str(max(delRange(:,2))) '.']);
            end
                
            % Save allocation parameters
            alloc.model = model;
            alloc.betaRange = betRange;
            alloc.deltaRange = delRange; 
            alloc.delOffset = delOffset;
            
            % Pre-compute cost function parameters 
            betaRangeInv = 1 ./ betRange;
            deltaRangeInv = 1 ./ (delOffset - delRange); 
            populationWeights = ones(alloc.model.n, 1);
            alloc.fSlope = populationWeights ./ (betaRangeInv(:,1) - betaRangeInv(:,2));
            alloc.gSlope = populationWeights ./ (deltaRangeInv(:,2) - deltaRangeInv(:,1)); 
            alloc.fOffset = alloc.fSlope .* betaRangeInv(:,2);
            alloc.gOffset = alloc.gSlope .* deltaRangeInv(:,1); 
        end
                
        function [R0, vaxAlloc, antiAlloc, newModel] = ...
                allocatePharmaR0(alloc, budget)
            % ALLOCATEPHARMAR0 Allocate vaccines and antidotes to minimize
            % R0. 
            % Inputs:
            %   budget      : upper bound on cost
            % Outputs:
            %   R0          : minimum value of R0
            %   vaxAlloc    : budget allocated to vaccines at each node
            %   antiAlloc   : budget allocated to antidotes at each node
            %   newModel    : post-intervention SEIR model
            
            n = alloc.model.n; 
            cvx_solver mosek
            cvx_precision default
            cvx_begin gp quiet

                variable r  
                variable w(2*n, 1) 
                variable bet(n, 1)
                variable del(n, 1)  % delta
                variable eta(n, 1)  % Proxy for delOffset - delta 

                % Construct NGM components
                F = [zeros(n) diag(bet.*alloc.model.s0)*alloc.model.A; ...
                    zeros(n, 2*n)];
                Vd = [diag(alloc.model.gamma) zeros(n); ...
                    zeros(n) diag(del)];
                Vod = [zeros(n, 2*n); diag(alloc.model.gamma) zeros(n)];

                % Construct total cost constraint
                f = sum(alloc.fSlope ./ bet); 
                g = sum(alloc.gSlope ./ eta);
                rhs = budget + sum(alloc.fOffset) + sum(alloc.gOffset);

                % Optimize 
                minimize(r)
                subject to   
                    diag(1/(r*Vd*w)) * (F+r*Vod) * w <= 1;  %#ok<*VUNUS>
                    f + g <= rhs;
                    alloc.betaRange(:,1) <= bet <= alloc.betaRange(:,2); 
                    alloc.delOffset - alloc.deltaRange(:,2) <= eta;
                    eta <= alloc.delOffset - alloc.deltaRange(:,1);
                    del + eta <= alloc.delOffset;   
            cvx_end
            
            % Unpack results 
            R0 = r;
            vaxAlloc = alloc.fSlope ./ bet - alloc.fOffset;
            antiAlloc = alloc.gSlope ./ eta - alloc.gOffset;
            newModel = SEIR(bet, alloc.model.gamma, ...
                alloc.delOffset-eta, alloc.model.s0, alloc.model.A);
            if ~strcmp(cvx_status, 'Solved')
                warning(['Solver returned status ' cvx_status]);
            end
        end 
        
        function [alpha, vaxAlloc, antiAlloc, newModel] = ...
                allocatePharmaAbscissa(alloc, budget)
            % ALLOCATEPHARMAABSCISSA Allocate vaccines and antidotes to
            % minimize spectral abscissa.
            % Inputs:
            %   budget      : upper bound on cost
            % Outputs:
            %   alpha       : minimum value of spectral abscissa
            %   vaxAlloc    : budget allocated to vaccines at each node
            %   antiAlloc   : budget allocated to antidotes at each node
            %   newModel    : post-intervention SEIR model 
            
            % Get diagonal offset, i.e., max_{\theta, i} (Vd)_i 
            % Offset ensures J + (diagOffset * eye(2*n)) is non-neg,
            % so we are guaranteed a positive eigenvalue
            J11max = max(alloc.model.gamma);
            J22max = max(alloc.deltaRange(:,2));
            if (alloc.delOffset < J11max) || (alloc.delOffset < J22max)
                error([...
                    'Offset ' num2str(alloc.delOffset) ...
                    ' is too small to transcribe the optimization' ...
                    ' as a geometric program. Try a value of' ...
                    ' delOffset that is at least ' ...
                    num2str(max(J11max, J22max))]);
            end 
             
            n = alloc.model.n; 
            cvx_solver mosek
            cvx_precision default
            cvx_begin gp quiet

                variable a              % -alpha
                variable w(2*n, 1) 
                variable bet(n, 1)      % beta  
                variable eta(n, 1)      % Proxy for delOffset - delta 
                
                % Construct non-neg shifted Jacobian matrix
                M11 = diag(alloc.delOffset - alloc.model.gamma);
                M12 = diag(bet.*alloc.model.s0)*alloc.model.A;
                M21 = diag(alloc.model.gamma);
                M22 = diag(eta);
                M = [[M11 M12]; [M21 M22]];

                % Construct total cost constraint
                f = sum(alloc.fSlope ./ bet); 
                g = sum(alloc.gSlope ./ eta);
                rhs = budget + sum(alloc.fOffset) + sum(alloc.gOffset);

                % Optimize
                minimize(a)
                subject to    
                    % Ensure a = abscissa(M)
                    M * w <= a * w;  
                    % Cost constraint
                    f + g <= rhs;
                    % Beta box constraint
                    alloc.betaRange(:,1) <= bet <= alloc.betaRange(:,2);
                    % Delta box constraint
                    alloc.delOffset - alloc.deltaRange(:,2) <= eta;
                    eta <= alloc.delOffset - alloc.deltaRange(:,1);
            cvx_end
            
            % Unpack results 
            alpha = a - alloc.delOffset;
            vaxAlloc = alloc.fSlope ./ bet - alloc.fOffset;
            antiAlloc = alloc.gSlope ./ eta - alloc.gOffset;
            newModel = SEIR(bet, alloc.model.gamma, ...
                alloc.delOffset - eta, alloc.model.s0, alloc.model.A);
            if ~strcmp(cvx_status, 'Solved')
                warning(['Solver returned status ' cvx_status]);
            end
        end
 
    end
end
