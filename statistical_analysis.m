class_data = readtable("../../Desktop/mcs_2021.xlsx");
%Selecting a specific column
 ode_col = class_data(:, {'ODE'}); 
 %disp(ode_col)
 %computing the mean value of each of the Units 
 ode_m = mean(class_data.("ODE"))
 computational_maths_m = mean(class_data.("ComputationalMath"))
real_analysis_m = mean(class_data.("RealAnalysis"))
fluid_mechanics_m = mean(class_data.("FluidMechanics"))
comp_graphics_m = mean(class_data.("ComputerGraphics"))
linear_algebra_m = mean(class_data.("LinearAlgebra"))
theory_of_est_m = mean(class_data.("TOE"))
numerical_analysis_m = mean(class_data.("NumericalAnalysis"))
%computing the mode of each of the units 
ode_md = mode(class_data.("ODE"))
computational_maths_md = mode(class_data.("ComputationalMath"))
real_analysis_md = mode(class_data.("RealAnalysis"))
fluid_mechanics_md = mode(class_data.("FluidMechanics"))
comp_graphics_md = mode(class_data.("ComputerGraphics"))
linear_algebra_md = mode(class_data.("LinearAlgebra"))
theory_of_est_md = mode(class_data.("TOE"))
numerical_analysis_md = mode(class_data.("NumericalAnalysis"))
%computing the median of each of the units 
ode_mdn = median(class_data.("ODE"))
computational_maths_mdn = median(class_data.("ComputationalMath"))
real_analysis_mdn = median(class_data.("REGNO_"))
fluid_mechanics_mdn = median(class_data.("FluidMechanics"))
comp_graphics_mdn = median(class_data.("ComputerGraphics"))
linear_algebra_mdn = median(class_data.("LinearAlgebra"))
theory_of_est_mdn = median(class_data.("TOE"))
numerical_analysis_mdn = median(class_data.("NumericalAnalysis"))
%calculating the mean mode and median for each student 
% lets first specify the column names to ensure we have everything
% correclty labelled
scoreColumns = {'ODE', 'ComputationalMath', 'RealAnalysis', 'FluidMechanics', 'ComputerGraphics', 'LinearAlgebra', 'TOE', 'NumericalAnalysis'};

% Preallocating the arrays to store the results
meanScores = zeros(height(class_data), 1);
medianScores = zeros(height(class_data), 1);
modeScores = cell(height(class_data), 1);

% Loop through each row and calculate mean, median, and mode for each student
for i = 1:height(class_data)
    scores = class_data{i, scoreColumns};
    meanScores(i) = mean(scores);
    medianScores(i) = median(scores);
    modeScores{i} = mode(scores);
end

% Display the results
results = table(class_data.SERIALNO_, class_data.REGNO_, class_data.NAMES, meanScores, medianScores, modeScores, 'VariableNames', {'SERIALNO_', 'REGNO_', 'NAMES', 'MeanScore', 'MedianScore', 'ModeScore'});
disp(results);

%Lets Test for uniformity and independence  
% Specify the column names
scoreColumns = {'ODE', 'ComputationalMath', 'RealAnalysis', 'FluidMechanics', 'ComputerGraphics', 'LinearAlgebra', 'TOE', 'NumericalAnalysis'};

% Preallocate arrays to store the test results
uniformityTestResults = zeros(length(scoreColumns), 2); % [p-value, reject H0]
independenceTestResults = zeros(length(scoreColumns) - 1, 2); % [p-value, reject H0]

% Looping through each score column and perform tests
for i = 1:length(scoreColumns)
    scores = class_data.(scoreColumns{i});
    
    % Performing Chi-square test for uniformity
    expectedUniform = ones(size(scores)) * mean(scores);
    [uniformityTestResults(i, 1), uniformityTestResults(i, 2)] = chi2gof(scores, 'Expected', expectedUniform);
    
    % Performing Chi-square test for independence
    if i < length(scoreColumns)
        nextScores = class_data.(scoreColumns{i+1});
        contingencyTable = crosstab(scores, nextScores);
        [independenceTestResults(i, 1), independenceTestResults(i, 2)] = chi2gof(contingencyTable(:));
    end
end

% Display the test results
uniformityResultsTable = table(scoreColumns', uniformityTestResults(:, 1), uniformityTestResults(:, 2), 'VariableNames', {'ScoreColumn', 'UniformityPValue', 'UniformityRejectH0'});
independenceResultsTable = table(scoreColumns(1:end-1)', independenceTestResults(:, 1), independenceTestResults(:, 2), 'VariableNames', {'ScoreColumn1', 'IndependencePValue', 'IndependenceRejectH0'});

disp('Uniformity Test Results:');
disp(uniformityResultsTable);

disp('Independence Test Results:');
disp(independenceResultsTable);

%Testing for correlation between each of the units 
% Specify the column names
scoreColumns = {'ODE', 'ComputationalMath', 'RealAnalysis', 'FluidMechanics', 'ComputerGraphics', 'LinearAlgebra', 'TOE', 'NumericalAnalysis'};

% Preallocate correlation matrix
correlationMatrix = zeros(length(scoreColumns));

% Calculate correlation coefficients
for i = 1:length(scoreColumns)
    for j = 1:length(scoreColumns)
        scores1 = class_data.(scoreColumns{i});
        scores2 = class_data.(scoreColumns{j});
        correlationMatrix(i, j) = corr(scores1, scores2);
    end
end

% Display correlation matrix
correlationResultsTable = array2table(correlationMatrix, 'VariableNames', scoreColumns, 'RowNames', scoreColumns);
disp('Correlation Matrix:');
disp(correlationResultsTable);




