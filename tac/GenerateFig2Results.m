
clear;
run('SetupExperiments');

% Initialize models
modelLow = SEIR(0.05 * ones(nNodes, 1), 0.2 * ones(nNodes, 1), ...
    0.2 * ones(nNodes, 1), s0, A);
modelMedium = SEIR(0.1 * ones(nNodes, 1), 0.2 * ones(nNodes, 1), ...
    0.1 * ones(nNodes, 1), s0, A);
modelHigh = SEIR(0.15 * ones(nNodes, 1), 0.2 * ones(nNodes, 1), ...
    0.075 * ones(nNodes, 1), s0, A);
allModels = {modelLow, modelMedium, modelHigh};

budgets = 0 : 0.05 : 10;
nBudgets = length(budgets);

disp(['Low model R0: ' num2str(modelLow.R0())]);
disp(['Medium model R0: ' num2str(modelMedium.R0())]);
disp(['High model R0: ' num2str(modelHigh.R0())]);
nModels = length(allModels);
 
% Prepare to store simulation / allocation results
baselineStats = newModelStats(nModels, nBudgets, nPts, nNodes);
minR0Stats = newModelStats(nModels, nBudgets, nPts, nNodes);
minAbscissaStats = newModelStats(nModels, nBudgets, nPts, nNodes);

% Create, simulate, and study all models
simNum = 1;
simCount = prod([nModels, nBudgets]);
for i = 1:nModels
    model = allModels{i};
    for j = 1:nBudgets
        budget = budgets(j);  

        disp(['Running simulation ' num2str(simNum) ' / ' num2str(simCount)]); 

        % Construct allocator
        allocBetaRange = [0.1 * model.beta model.beta];
        allocDeltaRange = [model.delta 2 * model.delta];
        alloc = SEIRAllocator(model, allocBetaRange, allocDeltaRange, delOffset); 

        % Write baseline stats 
        baselineStats = saveModelStats(baselineStats, i, j, model, x0, tf, tt);

        % Write R0-minimizing stats
        [R0, vaxAlloc, antiAlloc, modelR0] = ...
            alloc.allocatePharmaR0(budget);
        minR0Stats = saveModelStats(minR0Stats, i, j, modelR0, x0, tf, tt); 
        minR0Stats.vaxAlloc(i, j, :) = vaxAlloc;
        minR0Stats.antiAlloc(i, j, :) = antiAlloc; 
        
        % Write abscissa-minimizing stats
        [abscissa, vaxAlloc, antiAlloc, modelAbscissa] = ...
            alloc.allocatePharmaAbscissa(budget);
        minAbscissaStats = saveModelStats(minAbscissaStats, i, j, modelAbscissa, x0, tf, tt);
        minAbscissaStats.vaxAlloc(i, j, :) = vaxAlloc;
        minAbscissaStats.antiAlloc(i, j, :) = antiAlloc; 

        simNum = simNum + 1;
        
    end
end  

% Save results
save('./results/fixed-model-results');

function s = newModelStats(nModels, nBudgets, nPts, nNodes)
    s = struct;
    s.ifxTimeseries = zeros(nModels, nBudgets, nPts);
    s.R0 = zeros(nModels, nBudgets);
    s.Abscissa = zeros(nModels, nBudgets);
    s.CumIfx = zeros(nModels, nBudgets);
    s.PeakIfx = zeros(nModels, nBudgets);
    s.PeakTime = zeros(nModels, nBudgets);
    s.vaxAlloc = zeros(nModels, nBudgets, nNodes);
    s.antiAlloc = zeros(nModels, nBudgets, nNodes);
end

function statsStruct = saveModelStats(statsStruct, i, j, model, x0, tf, tt)

    % Run simulation
    [t, ~, e, z, r] = model.simulate(x0, tf); 

    % Timeseries of active infections
    % ode45 doesn't have fixed step sizes, so timeseries is 
    % interpolated to standard sample times tt
    ifxTotal = sum(e + z, 2);   
    ifxInterp = interp1(t, ifxTotal, tt, 'spline');  
    statsStruct.ifxTimeseries(i,j,:) = ifxInterp;
    
    % Basic stats
    statsStruct.R0(i,j) = model.R0();
    statsStruct.Abscissa(i,j) = model.abscissa();
    statsStruct.CumIfx(i,j) = ifxTotal(end) + sum(r(end, :));
    
    % Peak stats 
    [maxIfx, maxIdx] = max(ifxInterp);
    statsStruct.PeakIfx(i,j) = maxIfx;
    statsStruct.PeakTime(i,j) = tt(maxIdx); 
    
end
