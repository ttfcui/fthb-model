% This program computes the age distribution of first time home buyers. The
% model is explained in the accompanied PDF.

function [f_tau, mode, mean] = fthb_age_dist(p_H, theta, mu, sigma, save, bench)
    % p_H: House price
    % theta: Down payment requirement
    % mu: Drift on log wealth
    % sigma: Multiplier on BM

    a = log(theta*p_H)/sigma;
    c = mu/sigma;

    t = linspace(10^-4, 45, 1000);          % Time grid
    f_tau = a./((2*pi*t.^3).^(1/2)) .* ...  % Density
        exp(-(a-c*t).^2./(2*t));
    
    fig = figure;
    hold on
    if exist('bench', 'var')
        plot(20+t, bench, '--')
    end
    plot(20+t, f_tau, 'linewidth', 2);
    xlabel('Age', 'FontSize', 18)
    xlim(20 + [t(1), t(end)])
    ylabel('Density', 'FontSize', 18)
    set(gca, 'fontsize', 16)
    set(gcf, 'paperpositionmode', 'auto')
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4.5]; fig.PaperSize = [6 4.5];
    print('-dpdf', save)

    [~, max_ind] = max(f_tau);
    mode = 20+t(max_ind);

    mean = t*f_tau'/sum(f_tau) + 20;

    disp(['The mode is ' num2str(mode)])
    disp(['The mean is ' num2str(mean)])

end
