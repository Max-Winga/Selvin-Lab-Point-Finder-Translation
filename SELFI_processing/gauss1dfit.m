function [popt,ret] = gauss1dfit(ydata, xdata)

% Fitting 1d gaussian to data using lsqcurvefit in Matlab.
% There are two inputs. Both are required, to make it simple. The first is the raw
% data, the second is the x-axis data.
% An estimate will be calculated for initial parameters.
% The output contains the fitted data and the fitted parameters

format long;
% Unconstrained minimization
% fitting the gaussian model
% y(x)=p(1)+p(2)*exp(-(x-p(4))*(x-p(4))/(2*p(3)*p(3)) of gauss1d.m to noisy
% measurements

% Estimating the parameters. The center is obtained from the maximum.
% We are assuming that the peaks are all maxima and not minima. The
% standard deviation is obtained by finding the area under the ydata. It
% corresponds to an area that covers 68.3% of the total area from the center

% Need to substract off the constant term. Otherwise the total will be off,
% and so does the percentage. The easy way is to assume that the minimum is
% the constant term. Another way is perhaps take the average of the bottom
% 1% of the data. Or perhaps a better way is to take moving average of the
% whole gaussian, and then just take the smallest value. Let's do this.
% Smooth with 0.5% of the entire length of the data, with a minimum of 5.

% Finding the standard deviation
yLength = length(ydata);        % Finding the number of elements in ydata
if yLength < 6;                % Length of data must be more than 6
     error('gauss1dfit:lengthChk', 'Data must have more than 6 elements');
end

span = round(0.005*yLength);    % Finding the number that is 0.5% of yLength
if span < 5                     % If the span is less than 5, then make it 5
    span = 5;
end
ySmooth = smooth(ydata,span);   % Smoothing ydata with a span of 5 or more
[base,imin] = min(ySmooth);     % Finding the base of smoothed data by finding the minimum value

if imin-span>0 && imin+span<yLength     % Just making sure that matrix dimension is not exceeded
    base = mean(ySmooth((imin-span):(imin+span)));
end
yOffset = ydata - base;         % Substract ydata off the base to make the base zero

% keyboard
[ymax, imax] = max(yOffset);    % Find center
    if imax == 1 || imax == yLength
%         popt = zeros(1,7);
        popt = [];
        ret=[];
        return;
    end

% In order to find out the area under the curve, we need to take into
% account the width of each element of xdata. For that we can come out with
% a new array for xdata, that corresponds to the width of each point. The
% width of point 6, for example is simply the distance between point 7 and
% point 5, divided by 2. The width of the first and last point is simply
% the width of that point to the point next to it.

xdata1=[xdata(2)*2 xdata(3:end) xdata(end)*2];
xdata2=[xdata(1)*2 xdata(1:end-2) xdata(end-1)*2];
xwidth = (xdata1-xdata2)/2;

i = 1;                      % Initialize i from 1
yArea = yOffset.*xwidth;    % Get the area under the curve for each point
total = sum(yArea);         % Find the total summation of yOffset. The standard deviation covers 68.3% of total from center
StdevArea = 0.683*total;    % Find out the summation that is 68.3% from center
yStdev = yArea(imax);       % Initialize yStdev to the area under the peak

while (imax+i)<yLength && (imax-i)>1 && yStdev<StdevArea
    yStdev = yStdev + yArea(imax+i) + yArea(imax-i);
    i=i+1;
end
assignin('base','xdata',xdata);
assignin('base','ydata',ydata);
assignin('base','imax',imax);
% keyboard
Stdev = (xdata(imax+i)-xdata(imax-i))/2;

p0 = [base ymax Stdev xdata(imax)];  % Initial estimate: offset, amplitude, standard deviation and center
% keyboard
options = optimset('Display','off',...
    'Algorithm','levenberg-marquardt','Jacobian','on','TolFun',1e-6,'TolX',1e-6);
% keyboard
[popt,~,~,~,output] = lsqcurvefit(@gauss1dfunct,p0,xdata,ydata,[],[],options);
ret=output.iterations;
