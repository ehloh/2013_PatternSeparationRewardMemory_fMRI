function [cf_,gof] = f_eval2_createROCFit(XData,YData, graph)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(XDATA,YDATA)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1


% Data from dataset "ROC":
%    X = XData:
%    Y = YData:
%    Unweighted
%
% This function was automatically generated on 05-Oct-2010 11:50:25

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[580 102 680 485]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;


% --- Plot data originally in dataset "ROC"
XData = XData(:);
YData = YData(:);
h_ = line(XData,YData,'Parent',ax_,'Color',[0.333333 0 0.666667],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(XData));
xlim_(2) = max(xlim_(2),max(XData));
legh_(end+1) = h_;
legt_{end+1} = 'ROC';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
    xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
    set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[-0.0040000000000000000833, 0.40400000000000002576]);
end


% --- Create fit "Yonelinas"
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 0],'Upper',[3 1]);
ok_ = isfinite(XData) & isfinite(YData);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [0.5 0.2 ]; % D, R
set(fo_,'Startpoint',st_);
ft_ = fittype('x + R + ( (1-R) * normcdf( norminv(x,0,1) ,-D,1) ) - normcdf( norminv(x,0,1) ,0,1);',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'D', 'R'});

% Fit this model using new data
[cf_,gof] = fit(XData(ok_),YData(ok_),ft_,fo_);

% Or use coefficients from the original fit:
% if 0
%     cv_ = { 0.37628699836666446021, 0.34411950586135203745};
%     cf_ = cfit(ft_,cv_{:});
% end

% Plot this fit
h_ = plot(cf_,'predobs',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'Yonelinas Fit';
axis([0 1 0 1])

% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthWest'};
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
axis image
axis([0 1 0 1])
if graph==0
    close(f_)
end;
