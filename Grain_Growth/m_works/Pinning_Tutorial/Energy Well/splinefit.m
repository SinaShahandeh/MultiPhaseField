function splinefit(x,y)
%SPLINEFIT  Reconstruct figure in SPLINETOOL.
%
%   SPLINEFIT(X,Y) creates a plot, similar to the plot in SPLINETOOL,
%   using the data that you provide as input.
%   You can apply this function to the same data you used with
%   SPLINETOOL or with different data. You may want to edit the
%   function to customize the code or even this help message.

%   Make sure the data are in rows ...
x = x(:).'; y = y(:).';
% ... and start by plotting the data specific to the highlighted spline fit.

firstbox = [0.1300  0.4900  0.7750  0.4850];
subplot('Position',firstbox)
plot(x,y,'ok'), hold on
names={'data'};
ylabel('t_pos(:,2)')
xtick = get(gca,'Xtick');
set(gca,'xtick',[])

%   Now generate and plot the fit.

% spaps is used to enforce specified tolerance; :
spline1 = spaps(x,y,0.0028); 
names{end+1} = 'spline1'; fnplt(spline1,'-',2)


%   Plot the second graph from SPLINETOOL

subplot('Position',[ 0.1300  0.1300  0.7750  0.3100])
fnplt(fnder(spline1),2)
ylabel('1st Deriv. of spline1')
xlabel('t_pos(:,1)')


%   Return to plotting the first graph
subplot('Position', firstbox)


legend(names{:})
hold off
set(gcf,'NextPlot','replace')
