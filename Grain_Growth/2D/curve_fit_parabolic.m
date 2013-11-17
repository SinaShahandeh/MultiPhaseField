 function [n_fit, ci]=curve_fit_parabolic(timevec,r)
 
% --- Plot data originally in dataset "r vs. timevec"
timevec = timevec(:);
r = r(:);
plot(timevec,r,'.')

% --- Create fit "fit 2"
ok_ = isfinite(timevec) & isfinite(r);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
 
st_ = [1];
ft_ = fittype(['(k*t+' num2str(R0) '^(2))^0.5'],...
     'dependent',{'y'},'independent',{'t'},...
     'coefficients',{'k'});

% Fit this model using new data
%  cf_ = fit(timevec(ok_),r(ok_),ft_,'Startpoint',st_);
 cf_ = fit(timevec(ok_),r(ok_),ft_);
% Or use coefficients from the original fit:
if 0
   cv_ = { 1.4274637815643749228, 0.44102042725840889803, 14.662932663783200482};
   cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
% legend off;  % turn off legend from plot method call
% set(h_(1),'Color',[1 0 0],...
%      'LineStyle','-', 'LineWidth',2,...
%      'Marker','none', 'MarkerSize',6);
% legh_(end+1) = h_(1);
% legt_{end+1} = 'fit 2';

% Done plotting data and fits.  Now finish up loose ends.
% %hold off;
% leginfo_ = {'Orientation', 'vertical'}; 
% h_ = legend(ax_,legh_,legt_,leginfo_{:}); % create and reposition legend
% set(h_,'Units','normalized');
% t_ = get(h_,'Position');
% t_(1:2) = [0.633399,0.136249];
% set(h_,'Interpreter','none','Position',t_);
% xlabel(ax_,'');               % remove x label
% ylabel(ax_,'');               % remove y label
cf_
n_fit=cf_.n;
ci=confint(cf_,0.95);
ci=ci(:,2);
