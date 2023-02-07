function decayTerm = fitNegExp(X, Y)


% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
if size(Y,1) < size(Y,2)
    tbl = table(X', Y'); % row vectors after transpose (time,1)
else
    tbl = table(X',Y);
end

% Define the model as Y = a + exp(-b*x)
modelfun = @(b,x) b(1) + b(2) * exp(-b(3)*x(:, 1));  
beta0 = [0.5, 1, 1]; % Guess values to start with.  Just make your best guess.

try
mdl = fitnlm(tbl, modelfun, beta0);
catch
    'a'
end


coefficients = mdl.Coefficients{:, 'Estimate'};
yFitted = coefficients(1) + coefficients(2) * exp(-coefficients(3)*X);

try
decayTerm = coefficients(3);
catch
    'a'
end


% % plot data
% f = figure;
% plot(X, Y, 'b*', 'LineWidth', 2, 'MarkerSize', 8);
% grid on;
% % Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
% hold on;
% plot(X, yFitted, 'r-', 'LineWidth', 2);
% grid on;
% xlabel('X', 'FontSize', 10);
% ylabel('Y', 'FontSize', 10);
% legendHandle = legend('Noisy Y', 'Fitted Y', 'Location', 'north');
% legendHandle.FontSize = 25;


end
