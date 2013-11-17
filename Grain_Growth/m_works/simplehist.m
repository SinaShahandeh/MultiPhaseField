function [n,A]=simplehist(logGA)
%SIMPLEHIST    Create plot of datasets and fits
%   SIMPLEHIST(LOGGA)
%   Creates a plot, similar to the plot in the main distribution fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  0

% This function was automatically generated on 29-Nov-2008 16:21:55
 
% Data from dataset "logGA data":
%    Y = logGA
 
% Force all inputs to be column vectors
logGA = logGA(:);

% Set up figure to receive datasets and fits
% f_ = figure;
% figure(f_);
% set(f_,'Units','Pixels','Position',[441 134 680 465.975]);
% legh_ = []; legt_ = {};   % handles and text for legend
% ax_ = newplot;
% set(ax_,'Box','on');
% hold on;

% --- Plot data originally in dataset "logGA data"
t_ = ~isnan(logGA);
Data_ = logGA(t_);
[F_,X_] = ecdf(Data_,'Function','cdf'...
              );  % compute empirical cdf
Bin_.rule = 1;
[C_,E_] = dfswitchyard('dfhistbins',Data_,[],[],Bin_,F_,X_);

[N_,C_] = ecdfhist(F_,X_,'edges',E_); % empirical pdf from cdf

A=C_;
n=N_;
% h_ = bar(C_,N_,'hist');
% set(h_,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
%        'LineStyle','-', 'LineWidth',1);
% xlabel('Data');
% ylabel('Density')
% legh_(end+1) = h_;
% legt_{end+1} = 'logGA data';

% Nudge axis limits beyond data limits
% xlim_ = get(ax_,'XLim');
% if all(isfinite(xlim_))
%    xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
%    set(ax_,'XLim',xlim_)
% end

% x_ = linspace(xlim_(1),xlim_(2),100);

% hold off;
% leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
% h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
% set(h_,'Interpreter','none');
