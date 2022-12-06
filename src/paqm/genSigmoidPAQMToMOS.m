function [sigmoidPAQMToMOS, a1, a2, a3, a4] = genSigmoidPAQMToMOS(filepath)
    % GENERATESIGMOIDLNTOMOSCOEFFS computes the coefficients of the LnToMOS function.
    %
    % filepath must be some sort of table (like a .csv file).

    % Read values and convert table into vectors
    T = readtable(filepath);
    T = T{:, :};
    x = T(:, 1);
    y = T(:, 2);

    % Create the model for fitting
    modelfunc = @(b, x) b(3) + b(4) ./ (1 + exp(b(1) + b(2) .* x));
    options = statset('nlinfit');
    options.MaxIter = 10000;
    coeffs_vect = nlinfit(x, y, modelfunc, rand(1, 4), options);

    % Retrieve the values of the fitted model
    coeffs_cell = num2cell(coeffs_vect);
    [a1, a2, a3, a4] = coeffs_cell{:};

    sigmoidPAQMToMOS = @(p) a3 + a4 ./ (1 + exp(a1 + a2 * p));
    save sigmoidPAQMToMOS.mat a1 a2 a3 a4 sigmoidPAQMToMOS
end
