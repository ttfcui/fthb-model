% DESCRIPTION:
% This script loads the simulation results from fortran and plots the
% fthb steady state distribution and the transition dynamics...

% ARGUMENTS:
% POL: either -1 or -2. Describes if we want to extract information
% about first-time buyers (-1) or repeat buyers (-2).

% AUTHORS: GB, TC


% BEGIN
close all
% If calling as a script, can't let this delete input arguments!
% clc
% clear



policy_ss = importdata(strcat('dist_fthb','_',suffix,'.txt'));
hstock_orig = importdata(strcat('housingstock','_',suffix,'.txt'));
%Tretire = 38;                     % age of retirement (sample selection)
%End = 40;                         % number of periods
hstock_orig = hstock_orig(1:Tretire,:);
if POL == -1
    agType = 'FTHBs';
    fthb_ss = policy_ss((policy_ss(:,2) <= Tretire) & (policy_ss(:,17) == POL), :);
else
    agType = 'Repeat buyers';
    fthb_ss = policy_ss((policy_ss(:,2) <= Tretire) & (policy_ss(:,17) == POL), :);
end



%% Steady State Plots
% number of first time houses bought in every period

q_ss = size(fthb_ss(fthb_ss(:,2) <= Tretire, :), 1);
h_ss = sum(hstock_orig(:, 3)); % Total housing value transacted in steady state
t_ss = sum(hstock_orig(:, 4)) % Total housing transactions in steady state
i_ss = sum(fthb_ss(:, 3)); % Total *FTHB investment* in steady state

age = 21+ [1:Tretire];
fthb_age = zeros(Tretire, 1);
for t = 1:Tretire
    fthb_age(t) = sum(fthb_ss(:, 2)==t)/q_ss;
end

if POL == -1
% Age distribution (not a histogram because age is categorical)
% import actual FTHB age distribution
%try
    data_fthb = load('FTHB_dist_data.csv');
    f_frac = sprintf('%.1f', q_ss/t_ss*100);
    fig = figure;
        y = bar(age, [data_fthb(:,2), fthb_age], 'EdgeColor', 'w', 'BarWidth', 1);
        set(y(2), 'FaceColor', [1 0.6 0.6]);
        l = cell(1,2); l{1}='Data'; l{2}='Model';
        legend(y, l);
        set(gca, 'fontsize', 10)
        title('Age Distribution')
        ylabel(strjoin({'% of FTH bought by age (FTHBs ',f_frac,'% of all transactions)'}))
        xlabel('Age')
        xlim([21, age(end)]); ylim([0 0.08])
    
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4.5]; fig.PaperSize = [6 4.5];
    print('fracFthb','-dpdf')
    dlmwrite('FTHB_model_dist.txt',fthb_age,'delimiter','\t','precision',3)
%catch

% Size distribution
% Highlighting different income shock states
fig = figure;
    hold on
    colormap(fig, parula(4));
    twoway = scatter(fthb_ss(:, 5),fthb_ss(:, 3), 15, fthb_ss(:,8)./4, 'filled');
    twoway.MarkerEdgeAlpha = 0.0;
    twoway.MarkerFaceAlpha = 0.1;
    line = linspace(0, 3.0, 2); y = line;
    plot(line, y, '-r');
    betas = regress(fthb_ss(:, 3), [ones(size(fthb_ss, 1), 1), fthb_ss(:, 5)]);
    y = betas(1) + betas(2)*line; plot(line, y, '--k');
    xlabel('Rental housing size up to purchase')
    ylabel('House size upon purchase')
    xlim([0, 3.0]); ylim([0, 3.0]);
    % title('Histogram of fth size')
    set(gca, 'fontsize', 10)

fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4.5]; fig.PaperSize = [6 4.5];
% print('HouseTrans','-dpdf')

figure;
histogram(fthb_ss(:, 3), 50,'DisplayStyle','stairs','Normalization','probability');
xlabel('Size')
title('Histogram of fthb home sizes')
set(gca, 'fontsize', 16)
print('HouseSize', '-dpdf')

% Wealth distribution
figure;
histogram(fthb_ss(:, 7), 50,'DisplayStyle','stairs','Normalization','probability');
xlabel('Size')
title('Histogram of fthb wealth')
set(gca, 'fontsize', 16)
print('WealthHist', '-dpdf')

end

%% Transition to a policy change

% Looking only at working workers in transition
% (retirees have no income uncertainty so they exploit 1-period variation)
fthb_transition =  importdata(strcat('transition_fthb','_',suffix,'.txt'));
fthb_transition = fthb_transition((fthb_transition(:,3) <= Tretire) & ...
                                  (fthb_transition(:,18) == POL), :);
hstock = importdata(strcat('housing_transit','_',suffix,'.txt'));
hstock = hstock(hstock(:,2) <= Tretire, :);

PolYrs = 1; AppendStart = 0;
T_plot = End - 1;

% TODO: Append transition statistics when multi-period policies are in place

%% Affix all the period-cohort statistics together
Q_p = zeros(T_plot+1, 1); % FTHB proportion of buyers holder
Qf_p = zeros(T_plot+1, 1); % FTHB number holder
T_p = zeros(T_plot+1, 1); % Total transactions holder
H_p = zeros(End, 1); % Housing wealth holder
I_p = zeros(End, 1); % FTHB proportion of housing investment holder
If_p = zeros(End, 1); % FTHB housing investment holder

fthb_age_transition = zeros(Tretire, 2);
fthb_age_investment = zeros(Tretire, 2);
for p = 0:T_plot
    q_ss = size(fthb_ss(fthb_ss(:,2) > p, :), 1);
    p_ind = find((fthb_transition(:, 2) > 0) & ...
                 (fthb_transition(:, 3) - fthb_transition(:, 2) == p));
    h_ind = find(hstock(:,1) == p);
    fthb_p = fthb_transition(p_ind, :);
    Qf_p(p+1) = size(fthb_p,1)/q_ss;
    Q_p(p+1) = (size(fthb_p,1)/sum(hstock(h_ind,5))) - (q_ss/t_ss);
    if sum(fthb_p(:,8)) > 0
        p
        sum(fthb_p(:,8))/(Qf_p(p+1)*q_ss)
    end
    T_p(p+1) = sum(hstock(h_ind,5))/t_ss - 1;

    if (p < End)
    try
	%h_ind = find((hstock(:,1) == p) & (hstock(:,2) > T - p) & (hstock(:,2) <= T));
	owninvest_ss = hstock_orig(1:p, 2);
    catch % For period 0
        owninvest_ss = 0;
    end
    i_ss = sum(fthb_ss(fthb_ss(:,2) > p, 3));
    I_p(p+1) = sum(fthb_transition(p_ind, 4))/i_ss;
    If_p(p+1) = sum(fthb_transition(p_ind, 4));

    % TODO: Do we care about homeowners' housing wealth at all here?
    H_p(p+1) = (sum(fthb_transition(p_ind, 4)) + sum(owninvest_ss));
    end

end

% Percentage deviation in housing wealth, investment over last period
I_pd = I_p - (i_ss/h_ss);
If_pd =  If_p / i_ss;

% Get age distribution over time and investment ratios
q_ss = size(fthb_ss, 1);
for t = 1:Tretire
    p_ind = find(((fthb_transition(:, 3) - fthb_transition(:, 2)) < PolYrs) & ...
	 (fthb_transition(:, 3) == t));
    fthb_age_transition(t, 1) = size(fthb_transition(p_ind, 1), 1) ...
        /sum(Qf_p(1:PolYrs));
    % First column samples from policy-induced FTHBs,
    % Second samples from FTHBs in steady state
    fthb_age_investment(t, 1) = sum(fthb_transition(p_ind, 4))/ ...
        size(fthb_transition(p_ind, 1), 1);
    try
	fthb_age_investment(t, 2) = sum(fthb_ss(fthb_ss(:, 2) == t, 3))/ ...
	    size(fthb_ss(fthb_ss(:, 2) == t), 1);
    catch
        fthb_age_investment(t, 2) = 0; % Undefined due to no FTHBs of that age
    end

    p_ind = find(((fthb_transition(:, 3) - fthb_transition(:, 2)) >= PolYrs) & ...
	((fthb_transition(:, 3) - fthb_transition(:, 2)) < End) & ...
        (fthb_transition(:, 3) == t));
    fthb_age_transition(t, 2) = size(fthb_transition(p_ind, 1), 1) ...
        /sum(Qf_p(PolYrs+1:End));
end

for j = 1:2
    fthb_age_transition(1:Tretire, j) = ...
	fthb_age_transition(1:Tretire, j)/sum(fthb_age_transition(1:Tretire, j));
end

T_img = min(15,T_plot) ;  % reversal is done by this time
T_img
T_plot
fig = figure;
hold on
title('Transaction following a policy change')
T_pc = cumsum(T_p);
stairs(-3:T_img, [zeros(3, 1); T_p(1:T_img+1)])
stairs(-3:T_img, [zeros(3, 1); T_pc(1:T_img+1)])
plot(-3:T_img, zeros(T_img+4,1).', '-.k')
ylabel('Ratio of values relative to steady state transactions')
xlabel('Period')
% legend(strjoin({agType,'proportion among all transactions'}), ...
%        strjoin({agType,'housing investment prop. among all transactions'}),...
legend(strjoin({'Total transactions in period'}), ...
       strjoin({'Cum. transactions'}),...
       'Location','southoutside')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4.5]; fig.PaperSize = [6 4.5];
print('FthbShock','-dpdf')

lims = get(gca, 'YLim');
if (lims(1) > -0.25)
    lims(1) = -0.25;
end
if (lims(2) < 0.50)
    lims(2) = 0.50;
end

axis([get(gca, 'XLim'), lims])
print('FthbShockAxisAdj','-dpdf')



% Cumulative gains over steady-state flow
Q_pc = cumsum(Qf_p) - [1:T_plot+1]';
I_pc = I_p - 1;
fig = figure;
hold on
title('Transaction following a policy change')
stairs(-3:T_img, [zeros(3, 1); Q_pc(1:T_img+1)])
stairs(-3:T_img, [zeros(3, 1); I_pc(1:T_img+1)], 'r')
plot(-3:T_img, 0*(-3:T_img), '-.k')
ylabel(strjoin({'Ratio of values relative to steady state',agType}))
xlabel('Period')
legend(strjoin({'Cum.',agType,'purchases'}), ...
       strjoin({'Durable stocks, growth among',agType}),...
       'Location','southoutside')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4.5]; fig.PaperSize = [6 4.5];
print('FthbShockCum','-dpdf')
lims = get(gca, 'YLim');
if (lims(1) > -0.50)
    lims(1) = -0.50;
end
if (lims(2) < 1.00)
    lims(2) = 1.00;
end

axis([get(gca, 'XLim'), lims])
print('FthbShockCumAxisAdj','-dpdf')

dlmwrite('PolSeries.txt',[[zeros(3, 1); T_pc], [zeros(3, 1); Q_pc], ...
                          [zeros(3, 1); I_pc]],'delimiter',',','precision',6)

% Histograms of age distribution
fig = figure;
    subplot(2,1,1);
    hold on
    bar(age+0.5, fthb_age_transition(:,1))
    stairs(age, fthb_age, '-.', 'Color', 'r', 'LineWidth', 1.9)
    eval(['title(''Policy period( to period ' num2str(PolYrs) ')'');'])
    ylim([0 0.1])
    subplot(2,1,2);
    hold on
    bar(age+0.5, fthb_age_transition(:,2))
    stairs(age, fthb_age, '-.', 'Color', 'r', 'LineWidth', 1.9)
    title('Remaining Periods')
    xlabel('Age')
    ylabel(strjoin({'Prop. of',agType,'by age'}))
    ylim([0 0.1])
%     set(gca, 'fontsize', 16)
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4.5]; fig.PaperSize = [6 4.5];
print('FthbShockAge','-dpdf')
dlmwrite('FTHB_age_trans.txt',[[2:Tretire+1]', fthb_age_transition(:,1), fthb_age],...
         'delimiter','\t','precision',5)

% Histograms of policy period housing investment
fig = figure;
    hold on
    y = bar(age+0.1, fthb_age_investment, 'EdgeColor', 'w', 'BarWidth', 1);
    set(y(2), 'FaceColor', [1 0.6 0.6]);
    eval(['title(''Policy period( to period ' num2str(PolYrs) ')'');'])
    xlabel('Age')
    ylabel('Average housing size')
    ylim([0 6])
    xlim([20 60])
    l = cell(1,2); l{1}='Policy'; l{2}='Steady State';
    legend(y, l);
    set(fig, 'PaperOrientation','landscape', 'PaperUnits', 'inches', ...
        'PaperPosition', [1 1 10 6]);
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4.5]; fig.PaperSize = [6 4.5];
print('HouseInvShockAge','-dpdf')
dlmwrite('AgeInvLevels.txt',fthb_age_investment,'delimiter','\t','precision',3)

% plot(age, fthb_age, '--', age, fthb_age_transition(:, 1), age, fthb_age_transition(:, 5), 'linewidth', 3)
% legend('Steady State', 'On impact', 'After 4 periods')

% END
