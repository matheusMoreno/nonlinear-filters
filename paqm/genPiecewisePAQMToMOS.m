function [piecewisePAQMToMOS, a, b, c, d, beta_1, alpha_1, beta_2, alpha_2] = genPiecewisePAQMToMOS(filepath, p_1, p_2)
    % GENPIECEWISEPAQMTOMOS generates the piecewise LnToMOS function.

    % Retrieve data points
    table = readtable(filepath);
    table = table{:, :};
    x = table(:, 1);
    y = table(:, 2);

    % Fit a third-order polynomial
    coefs = polyfit(x, y, 3);

    % Retrieve the values of the fitted model
    coefs_cell = num2cell(coefs);
    [a, b, c, d] = coefs_cell{:};

    assert(                                                         ...
        (-2 * b - sqrt(4 * b ^ 2 - 12 * a * c)) / (6 * a) <= p_1,   ...
        'p_1 is smaller than the local maximum.'                    ...
    );

    assert(                                                         ...
        (-2 * b + sqrt(4 * b ^ 2 - 12 * a * c)) / (6 * a) >= p_2,   ...
        'p_2 is bigger than the local minimum.'                     ...
    );

    % Compute terms for f_1(x)
    y_1 = polyval(coefs, p_1);
    beta_1 = abs(3 * a * p_1 .^ 2 + 2 * b * p_1 + c) ./ (5 - y_1);
    alpha_1 = (5 - y_1) ./ (exp(p_1 * beta_1));

    % Compute terms for f_3(x)
    y_2 = polyval(coefs, p_2);
    beta_2 = abs(3 * a * p_2 .^ 2 + 2 * b * p_2 + c) ./ (y_2 - 1);
    alpha_2 = (y_2 - 1) ./ (exp(-p_2 * beta_2));

    % fprintf('a = %.4f, b = %.4f, c = %.4f, d = %.4f\n', a, b, c, d);
    % fprintf('beta_1 = %.4f, alpha_1 = %.4f, beta_2 = %.4f, alpha_2 = %.4f\n', ...
    %         beta_1, alpha_1, beta_2, alpha_2);

    function m = pwfunc(p)
        if p <= p_1
            m = 5 - alpha_1 * exp(beta_1 .* p);     % f_1
        elseif p <= p_2
            m = polyval([a b c d], p);              % f_2
        else
            m = 1 + alpha_2 * exp(-beta_2 .* p);    % f_3
        end
    end

    piecewisePAQMToMOS = @(p) pwfunc(p);
    save piecewisePAQMToMOS.mat a b c d beta_1 alpha_1 beta_2 alpha_2 p_1 p_2 piecewisePAQMToMOS
end
