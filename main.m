function main(FNAME, ALPHAS, ISMALE)

%MAIN analysis file to accompany the paper:
%
%  ANGUS, SD: 'A statistical timetable for the sub-2 hour marathon'
%  Dept. of Economics, Monash University
%  to appear in MSSE
%
%  MAIN(FNAME,ALPHAS,ISMALE) calculates all necessary measures
%     and produces figures, using input CSV file FNAME, for
%     chosen ALPHAS vector, and a switch for whether or not
%     the FNAME provides male (ISMALE=1) or female(ISMALE=0)
%     data (this changes the pathway for calculation and figure
%     generation).
%
%  OUTPUTS
%     When run, MAIN will produce two gendergap files for later
%     comparison with gender_gap(). For instance, if ISMALE=1
%     the two files will be: out_expected_male.csv, and
%     out_1in10_male.csv . After running both male and female
%     variants, one can then conduct a long-run gender gap analysis
%     as follows:
%     >> OUTNAME = 'gg_1in10.csv';
%     >> gender_gap('out_1in10_male.csv','out_1in10_female.csv',OUTNAME)
%
%  Notes:
%   -- tested with MATLAB R2018b
%   -- requires Statistics Toolbox installed
%
%   See also FITNLM GEVFIT PREDICT

SUB2_MIN = 120;

SUB2_CROSSING_POINT = SUB2_MIN - 1/60;  % less 1s to ensure the barrier is broken
ALPHAS_INSET = linspace(0.005,0.5)';    % alpha values for inset figure
ODDS_FOCUS = 10;                        % benchmark 'odds' (10% likely) for baseline analysis

% .. check license
if license('checkout', 'statistics_toolbox') == 0
    error('Did not Statistics Toolbox installed. Exiting.')
end

% // PREP DATA ------------------- %
if ISMALE, SUF = 'male', else SUF = 'female', end;  % .. for output file endings
SUB2_CP_pc = SUB2_CROSSING_POINT / SUB2_MIN;        % percentage terms for -1s
def_clrs = get(groot,'defaultAxesColorOrder');      % store matlab defaults

% .. ingest
db = readtable(FNAME, 'delimiter', '|');

% .. drop unofficial
db = db(db.unofficial ~= 1,:);

% .. drop years before 1950
db = db(db.date.Year >= 1950,:);

fprintf(' --> after filtering, we have %.0f observations.\n', height(db))

% .. obtain useful vars
db.dt  = datetime(db.time);
db.hrs = db.dt.Hour;
db.min = db.dt.Minute;
db.sec = db.dt.Second;
db.tot_min = db.hrs*60 + db.min + db.sec/60;

% .. index relative SUB2_MIN setting (=1.00)
db.sub2 = db.tot_min/SUB2_MIN;

% .. get years since first date as X
db.date_dt = datetime(db.date);
date0 = db.date_dt(1);
db.date_diff_yrs = years(db.date_dt - date0);   % a duration

% // MODEL ------------------- %
X = db.date_diff_yrs;
Y = db.sub2;
fun = @(beta,x)(beta(1)*exp(-beta(2)*x) + beta(3));     % T = B0.e^(-lambda.y) + B_inf
beta0 = [1 0.1 1.1];      % initial conditions
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
nlm = fitnlm(X,Y,fun,beta0, 'Options', opts, 'CoefficientNames', {'B0' 'lambda' 'B_inf'});
% .. display
disp(nlm)

% // extension: fit residuals to GEV distribution, per Denny2008
fprintf(' --> GEV fit to -ve of residuals ... \n')
opts = statset('MaxIter', 1E5);
x = -nlm.Residuals.Raw;     % neg, as we want max lower residual value
alpha = 0.05;
[params,ci] = gevfit(x,alpha,opts);
disp(params)
disp(ci)
if params(1) < 0
    maxV = params(2) - params(3)/params(1);
    fprintf(' --> a < 0, maxV ~ %.5f\n',maxV)
else
    fprintf(' --> a >= 0, no max specified.\n')
end


% // Predictions ------------------- %
% .. for main plot
x_pred = linspace(-20,150)';
for i = 1:numel(ALPHAS)
    pred(i).ALPHA = ALPHAS(i);
    [pred(i).y,pred(i).ci] = predict(nlm, x_pred, 'Prediction', 'observation', 'Alpha', pred(i).ALPHA, 'Simultaneous', true);
end

% .. for inset
if ISMALE
    for i = 1:numel(ALPHAS_INSET)
        pred_inset(i).ALPHA = ALPHAS_INSET(i);
        [~,pred_inset(i).ci] = predict(nlm, x_pred, 'Prediction', 'observation', 'Alpha', pred_inset(i).ALPHA, 'Simultaneous', true);
        t = interp1(pred_inset(i).ci(:,1), x_pred, SUB2_CP_pc);
        pred_inset(i).date_sub2 = date0 + t*365;
    end
    a1 = 1./([pred_inset.ALPHA]/2);
    p1 = dt2year_frac([pred_inset.date_sub2]);
end


% // ANALYSIS ------------------- %
% Q1: WHEN WILL A RUNNER GO SUB-2?
% .. interp date at which obs. prediction will == SUB2_CROSSING_POINT;
if ISMALE
    fprintf('\nSUB2\n')
    for i = 1:numel(ALPHAS)
        t = interp1(pred(i).ci(:,1), x_pred, SUB2_CP_pc);
        if ~isnan(t)
            pred(i).date_sub2 = date0 + t*365;
            fprintf(' --> sub2 date (alpha=%4.3f, 1 in %.0f): %s\n', pred(i).ALPHA, 1/(pred(i).ALPHA/2), pred(i).date_sub2);
        else
            fprintf(' --> sub2 date (alpha=%4.3f): - no crossing -\n', pred(i).ALPHA);
        end
    end
end

% Q2: WHAT IS THE LONG-RUN BEST PERFORMANCE POSSIBLE?
% .. obtain run-time when year-->\inf, compute for various levels of odds
for i = 1:numel(ALPHAS_INSET)
    pred_best(i).ALPHA = ALPHAS_INSET(i);
    x_best = inf;       % infinite year
    [~,pred_best(i).ci] = predict(nlm, x_best, 'Prediction', 'observation', 'Alpha', pred_best(i).ALPHA, 'Simultaneous', true);
    pred_best(i).limit_time = pred_best(i).ci(1) * SUB2_MIN;    % convert back to min
end
a2 = 1./([pred_best.ALPHA]/2);
p2 = [pred_best.limit_time];

% Q3: WHAT IS THE CURRENT PERFORMANCE GAP (FROM WR TO BEST POSSIBLE)?
% .. obtain differnce between current_wr and limit, express as % difference
fprintf('\nPERFORMANCE GAP\n')
curr_wr = db.tot_min(end);
for i = 1:numel(ALPHAS)
    gap(i).alpha = ALPHAS(i);
    gap(i).odds  = 1 ./ (ALPHAS(i)/2);
    gap(i).best_min = interp1(a2,p2, gap(i).odds);
    gap(i).gap_min  = round(gap(i).best_min - curr_wr, 2);
    gap(i).gap_pc   = round(gap(i).gap_min./curr_wr * 100,2);
end
table_gap = struct2table(gap);
disp(table_gap)

% // FIGURES ------------------- %

% figure 1
figure(1),clf, set(gcf,'Color','w')
    set(gcf,'Position',[1     1   549   672])
    hold on
    x_viz = db.date_diff_yrs*365 + date0;
    % .. main curve
    x = x_pred*365 + date0;
    y = round(pred(1).y,4);
    plot(x, y, '-', 'LineWidth', 3, 'Color', [0 0 0 0.5])
    % .. // export for gender gap analysis
    writetable(table(x,y),sprintf('out_expected_%s.csv',SUF));
    % .. add model predictions
    for i = 1:numel(ALPHAS)
        % .. main interest: lower tail
        plot(x_pred*365 + date0, pred(i).ci(:,1), 'r-', 'LineWidth', 1)
        % .. upper tail
        plot(x_pred*365 + date0, pred(i).ci(:,2), '-', 'LineWidth', 0.5, 'Color', [0 0 0 0.5])
        % .. sub2 points if ISMALE
        if ISMALE & ~isempty(pred(i).date_sub2)
            plot(pred(i).date_sub2, SUB2_CP_pc, 'k.', 'MarkerSize', 12)
        end
        % .. // export data for gender gap analysis
        if pred(i).ALPHA==0.2       % 1 in 10
            y = round(pred(i).ci(:,1),4);
            writetable(table(x,y),sprintf('out_1in10_%s.csv',SUF));
        end
    end
    if ~ISMALE      % add 'sub-130' line
        xx = [date0 date0+1*365];
        plot(xx, 130/SUB2_MIN * [1 1],'k--','LineWidth',0.5)
    end
    % .. add points last
    scatter(x_viz, db.sub2, 40, def_clrs(1,:), 'filled')
    % .. pretty
    if ISMALE
        set(gca,'Xscale', 'linear', 'Xlim', [datetime(1945,0,0) datetime(2105,0,0)])
    else
        set(gca,'Xscale', 'linear', 'Xlim', [datetime(1960,0,0) datetime(2020,0,0)])
    end
    set(gca,'FontSize',16)
    xlabel('Year')
    ylabel('Elapsed Time (h:min)')
    grid on
    % .. convert back to minutes on y-axis
    y_ticks = get(gca,'YTick');
    set(gca,'YTickLabel', nice_ticks(y_ticks*SUB2_MIN))
    % .. save/print
    fname = strrep(FNAME, '.csv', '_main-fig.eps');
    print(gcf,fname, '-depsc')

% figure 1: inset
if ISMALE
figure(2),clf,set(gcf,'Color','w')
    set(gcf,'Position', [244   416   307   257])
    plot(a1,p1,'-', 'LineWidth', 2), hold on
    box off, grid on
    set(gca,'XScale', 'log')
    set(gca, 'FontSize', 14)
    % .. get crossing points for annotation
    fprintf('\nINSET: ANNOTATION\n')
    inset_years = [2010 2020 2030 2040 2050];
    for i = 1:numel(inset_years)
        odds(i,1) = interp1(p1,a1,inset_years(i));
       fprintf(' --> inset odds at year %.0f: 1 in %.0f\n', inset_years(i), odds(i))
    end
    scatter(odds, inset_years', 40, 'filled')
    % .. get point at ODDS_FOCUS
    odds_focus_year = interp1(a1,p1,ODDS_FOCUS);
    fprintf(' --> inset year at odds-focus (1 in %.0f): %.0f\n', ODDS_FOCUS, odds_focus_year)
    plot(ODDS_FOCUS, odds_focus_year, 'k+', 'MarkerSize', 14)
    % ..dress
    xlabel('Odds of sub-2')
    ylabel('Year')
    % .. convert back to minutes on y-axis
    y_ticks = get(gca,'YTick');
    set(gca,'YTickLabel', nice_ticks(y_ticks))
    fname = strrep(FNAME, '.csv', '_sub2-inset.eps');
    print(gcf,fname, '-depsc')
end

% figure 2: best possible
figure(3),clf,set(gcf,'Color','w')
    set(gcf,'Position', [1   382   368   291])
    plot(a2,p2,'-', 'LineWidth', 2), hold on
    box off, grid on
    set(gca,'XScale', 'log')
    set(gca, 'FontSize', 14)
    % .. get crossing points for annotation
    fprintf('\nBEST: ANNOTATION\n')
    if ISMALE
        best_times = [116:119];
    else
        best_times = [122:2:128];
    end
    for i = 1:numel(best_times)
        odds_best(i,1) = interp1(p2,a2,best_times(i));
       fprintf(' --> best odds at time %.0fmin: 1 in %.0f\n', best_times(i), odds_best(i))
    end
    scatter(odds_best, best_times', 40, 'filled')
    % .. get point at ODDS_FOCUS
    odds_focus_limit = interp1(a2,p2,ODDS_FOCUS);
    fprintf(' --> limit time at odds-focus (1 in %.0f): %.2fmin\n', ODDS_FOCUS, odds_focus_limit)
    plot(ODDS_FOCUS, odds_focus_limit, 'k+', 'MarkerSize', 14)

    % .. dress
    xlabel('Odds of Limiting Time')
    ylabel('Limiting Time')
    y_ticks = get(gca,'YTick');
    set(gca,'YTickLabel', nice_ticks(y_ticks))
    fname = strrep(FNAME, '.csv', '_limit-time.eps');
    print(gcf,fname, '-depsc')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yf = dt2year_frac(dt)
yf = dt.Year + dt.Month/12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = nice_ticks(v)
% v: in minutes

for i = 1:numel(v)
    tot_m = v(i);
    h = floor(tot_m/60);
    m = num2str(floor(tot_m - 60*h));
    % .. ensure minutes has leading zero where req.
    if numel(m)<2
        m = ['0' m];
    end
    c{i} = sprintf('%.0f:%s',h,m);
end
