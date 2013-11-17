function [A1,A2,Delu]=ExtractDist(A,n,GrainsA,doplot)
%% Extract two size population density Based on number density
% A curve fitting code, fits double Gaussian fucntion to distribution
if strcmp(doplot,'yes')
    f_ = gcf;
%     figure(f_);
    set(f_,'Units','Pixels','Position',[441 134 680 475]);
    legh_ = []; legt_ = {};   % handles and text for legend
    xlim_ = [Inf -Inf];       % limits of x axis
    ax_ = axes;
    set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
    set(ax_,'Box','on');
    axes(ax_); hold on;
end
% --- Plot data originally in dataset "y vs. x"
x = A(:);
y = n(:);
% y=n(:)/sum(n(:));
xlim_(1) = min(x)-0.2*min(x);
xlim_(2) = max(x)+0.1*max(x);
if strcmp(doplot,'yes')
    h_ = line(x,y,'Parent',ax_,'Color',[0.333333 0 0.666667],...
        'LineStyle','none', 'LineWidth',1,...
        'Marker','.', 'MarkerSize',12);
    legh_(end+1) = h_;
    legt_{end+1} = 'y vs. x';
    % Nudge axis limits beyond data limits
    if all(isfinite(xlim_))
        xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
        set(ax_,'XLim',xlim_)
    else
        set(ax_, 'XLim',[-15.161660148229525902, -3.0549015239848000824]);
    end
end

% --- Create fit "fit 2"
ave_=mean(log(GrainsA));
std_=std(log(GrainsA));

coeffL=[0 0 0 0 ave_ 0];
coeffU=[1 ave_ std_ 1 ave_+2*std_ std_]; 
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',...
    coeffL,'Upper',coeffU);
ok_ = isfinite(x) & isfinite(y);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [0.1 ave_-1 1 0.1 ave_+1 1 ];
set(fo_,'Startpoint',st_);
ft_ = fittype('gauss2');
% Fit this model using new data
cf_ = fit(x(ok_),y(ok_),ft_,fo_);

coeff=coeffvalues(cf_);
xi=linspace(xlim_(1),xlim_(2),1000);
y1=coeff(1)*exp(-((xi-coeff(2))/coeff(3)).^2);
y2=coeff(4)*exp(-((xi-coeff(5))/coeff(6)).^2);

if strcmp(doplot,'yes')
    % Plot this fit
    h_ = plot(cf_,'fit',0.95);
    legend off;  % turn off legend from plot method call
    set(h_(1),'Color',[1 0 0],...
        'LineStyle','-', 'LineWidth',2,...
        'Marker','none', 'MarkerSize',6);
    legh_(end+1) = h_(1);
    legt_{end+1} = 'Gaussian Fit';
    % Done plotting data and fits.  Now finish up loose ends.
    hold off;
    leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
    h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
    set(h_,'Interpreter','none');
    xlabel(ax_,'');               % remove x label
    ylabel(ax_,'');               % remove y label
    % plotting individual distributions
    hold on

    plot(xi,y1,'g-');
    plot(xi,y2,'b-.');
    grid on
end
%% Analysis of volume fraction
A1=trapz(xi,y1);
A2=trapz(xi,y2);

Delu=coeff(5)-coeff(2);
hold off
pause(0.01);

