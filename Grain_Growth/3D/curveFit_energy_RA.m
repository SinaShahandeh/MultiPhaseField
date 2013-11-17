function curveFit_energy_RA(finalsize,TotalE)
%CURVEFIT_ENERGY_RA    Create plot of datasets and fits
%   CURVEFIT_ENERGY_RA(FINALSIZE,TOTALE)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1

 
% Data from dataset "TotalE vs. finalsize":
%    X = finalsize:
%    Y = TotalE:
%    Unweighted
%
% This function was automatically generated on 17-Aug-2011 18:09:38

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[998 258 674 477]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;

 
% --- Plot data originally in dataset "TotalE vs. finalsize"
finalsize = finalsize(:);
TotalE = TotalE(:);
h_ = line(finalsize,TotalE,'Parent',ax_,'Color',[0.333333 0 0.666667],...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(finalsize));
xlim_(2) = max(xlim_(2),max(finalsize));
legh_(end+1) = h_;
legt_{end+1} = 'TotalE vs. finalsize';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[33.536468448181629753, 121.28378424546319536]);
end


% --- Create fit "fit 1"
ok_ = isfinite(finalsize) & isfinite(TotalE);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [0.62274381989476468302 ];
ft_ = fittype('300^3*em-4/3*pi*em*x^3+4*pi*x^2*0.6667/2',...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'em'});

% Fit this model using new data
cf_ = fit(finalsize(ok_),TotalE(ok_),ft_,'Startpoint',st_);

% Or use coefficients from the original fit:
if 0
   cv_ = { 0.090460690284568531783};
   cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
     'LineStyle','-', 'LineWidth',2,...
     'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 1';

% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
