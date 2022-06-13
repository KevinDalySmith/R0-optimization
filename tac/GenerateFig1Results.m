
clear;
run('SetupExperiments');

% Specify ranges for baselines 
betaVals = 0.025 : 0.025 : 0.5;
gammaVals = 0.05 : 0.05 : 0.5;
deltaVals = 0.05 : 0.05 : 0.5;  
[betas, gammas, deltas] = meshgrid(betaVals, gammaVals, deltaVals); 

% Prepare to store simulation / allocation results
nBeta = length(betaVals);
nGamma = length(gammaVals);
nDelta = length(deltaVals);
baselineStats = newModelStats(nBeta, nGamma, nDelta, nPts, nNodes);
minR0Stats = newModelStats(nBeta, nGamma, nDelta, nPts, nNodes);
minAbscissaStats = newModelStats(nBeta, nGamma, nDelta, nPts, nNodes);

% Create, simulate, and study all models
simNum = 1;
simCount = prod([nBeta, nGamma, nDelta]);
for i = 1:nBeta
    beta = betaVals(i) * ones(nNodes, 1); 
    for j = 1:nGamma
        gamma = gammaVals(j) * ones(nNodes, 1);
        for k = 1:nDelta
            delta = deltaVals(k) * ones(nNodes, 1);
            
            disp(['Running simulation ' num2str(simNum) ' / ' num2str(simCount)]);

            % Construct model
            model = SEIR(beta, gamma, delta, s0, A);   
            
            % Construct allocator
            allocBetaRange = [0.1 * beta beta];
            allocDeltaRange = [delta 2 * delta];
            alloc = SEIRAllocator(model, allocBetaRange, allocDeltaRange, delOffset); 
            
            % Write baseline stats 
            baselineStats = saveModelStats(baselineStats, i, j, k, model, x0, tf, tt);
            
            % Write R0-minimizing stats
            [R0, vaxAlloc, antiAlloc, modelR0] = ...
                alloc.allocatePharmaR0(budget);
            minR0Stats = saveModelStats(minR0Stats, i, j, k, modelR0, x0, tf, tt); 
            minR0Stats.vaxAlloc(i, j, k, :) = vaxAlloc;
            minR0Stats.antiAlloc(i, j, k, :) = antiAlloc; 
            
            % Write abscissa-minimizing stats
            [abscissa, vaxAlloc, antiAlloc, modelAbscissa] = ...
                alloc.allocatePharmaAbscissa(budget);
            minAbscissaStats = saveModelStats(minAbscissaStats, i, j, k, modelAbscissa, x0, tf, tt);
            minAbscissaStats.vaxAlloc(i, j, k, :) = vaxAlloc;
            minAbscissaStats.antiAlloc(i, j, k, :) = antiAlloc; 
            
            simNum = simNum + 1;
 
        end
    end
end  

% Save results
save('./results/fixed-budget-results');

function s = newModelStats(nBeta, nGamma, nDelta, nPts, nNodes)
    s = struct;
    s.ifxTimeseries = zeros(nBeta, nGamma, nDelta, nPts);
    s.R0 = zeros(nBeta, nGamma, nDelta);
    s.Abscissa = zeros(nBeta, nGamma, nDelta);
    s.CumIfx = zeros(nBeta, nGamma, nDelta);
    s.PeakIfx = zeros(nBeta, nGamma, nDelta);
    s.PeakTime = zeros(nBeta, nGamma, nDelta); 
    s.vaxAlloc = zeros(nBeta, nGamma, nDelta, nNodes);
    s.antiAlloc = zeros(nBeta, nGamma, nDelta, nNodes);
end

function statsStruct = saveModelStats(statsStruct, i, j, k, model, x0, tf, tt)

    % Run simulation
    [t, ~, e, z, r] = model.simulate(x0, tf); 

    % Timeseries of active infections
    % ode45 doesn't have fixed step sizes, so timeseries is 
    % interpolated to standard sample times tt
    ifxTotal = sum(e + z, 2);   
    ifxInterp = interp1(t, ifxTotal, tt, 'spline');  
    statsStruct.ifxTimeseries(i,j,k,:) = ifxInterp;
    
    % Basic stats
    statsStruct.R0(i,j,k) = model.R0();
    statsStruct.Abscissa(i,j,k) = model.abscissa();
    statsStruct.CumIfx(i,j,k) = ifxTotal(end) + sum(r(end, :));
    
    % Peak stats 
    [maxIfx, maxIdx] = max(ifxInterp);
    statsStruct.PeakIfx(i,j,k) = maxIfx;
    statsStruct.PeakTime(i,j,k) = tt(maxIdx); 
    
end
