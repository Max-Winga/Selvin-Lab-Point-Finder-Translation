function [popt,resnorm,residual,ret]=gauss2dfit(zdata,grid,tilt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fitting 2d gaussian function onto a 2d array of data ('zdata') using
    % the lsqcurvefit (Levenberg-Marquadt) algorithm in Matlab
    
    % Input: 
    % 'zdata' is the raw two dimensional data
    % 'grid' is an array consisting of the dimension in x (grid(1)) and y (grid(2)),
    % and the gridding in x and y (which will be placed right after the 
    % dimensions, in the same column array). This input is [nx ny x y].
    % 'tilt' is an optional input. If the value is one, fitting will be
    % done by considering the tilt of the 2d gaussian

    % Output:
    % The output consists of the fitted parameters ('popt'), the number of
    % iterations ('ret') and the initial parameters used ('pInit')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Example
    % data = ones(100,100);               % Build data
    % dimension=size(data);               % Find out the dimensions
    % [ny,nx]=size(data);                 % Find out the dimensions
    % xAxis = 1:nx;                       % Create the x-Axis gridding
    % yAxis = 1:ny;                       % Create the y-Axis gridding
    % grid = [nx ny xAxis yAxis];         % Gridding input for gauss2dfunct and gauss2dfit
    % zOffset = 4.0;                      % Create parameters
    % Amplitude = 100.0;
    % xStd = 10.0;
    % yStd = 15.0;
    % xCenter = 45;
    % yCenter = 55;
    % tilt=0;                             % Value: 0 or 1. Specify whether or not there is tilt
    % tiltVal=0;                          % Value of tilt if tilt is to be evaluated
    % p0=[zOffset, Amplitude, xStd, yStd, xCenter, yCenter,tiltVal];      % Combine parameters
    % z = gauss2dfunct(p0,grid);            % Find the gaussian function
    % z = reshape(z,ny,nx);             % Convert 1d array into 2d array
    % z = z+Amplitude/10*rand(ny,nx);   % Add some random noise
    % [popt,ret]=gauss2dfit(z,grid,tilt);     % Run gauss2dfit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %tic;

    % Determine the x and y dimensions of data
    nx=grid(1);
    ny=grid(2);
    
    % Determine the dimensions of zdata
    dim=size(zdata);
    
    % Check that the dimensions given agree with dimensions of zdata
    if nx ~= dim(2)
        error('gauss2dfit:dimChk', 'x-dimension given must agree with x-dimension of raw data');
    end
    if ny ~= dim(1)
        error('gauss2dfit:dimChk', 'y-dimension given must agree with y-dimension of raw data');
    end

    % Convert 2d array into 1d array
    z = reshape(zdata,ny*nx,1);  
%     keyboard
    % Determine the peak
    [~, imax] = max(z);
    
    %ix = mod(imax,nx)+1;          % Find the x element of peak
    %iy = floor(imax/nx)+1;        % Find the y element of peak
    %xAxis = grid(3:nx+2);     % Get the x-axis for 1d gaussian fitting
    %yAxis = grid(nx+3:end);   % Get the y-axis for 1d gaussian fitting 
    
    %[px, ~] = gauss1dfit(zdata(ix,:),xAxis);  % 1d gaussian fit to the x-axis of data
    %[py, ~] = gauss1dfit(zdata(:,iy)',yAxis);  % 1d gaussian fit to the y-axis of data
%     keyboard
    ix = floor(imax/ny)+1;          % Find the x element of peak
    iy = mod((imax-1),ny)+1;        % Find the y element of peak

    if ix > nx
        ix = nx;
    end

    if iy > ny
        iy = ny;
    end
    
    xAxis = grid(3:nx+2);     % Get the x-axis for 1d gaussian fitting
    yAxis = grid(nx+3:end);   % Get the y-axis for 1d gaussian fitting 
    
    assignin('base','zdata',zdata);
%     keyboard
    [px, ~] = gauss1dfit(zdata(iy,:),xAxis);  % 1d gaussian fit to the x-axis of data
    [py, ~] = gauss1dfit(zdata(:,ix)',yAxis);  % 1d gaussian fit to the y-axis of data    
    
    % Find the initial estimate of the parameters
    % For 1d gaussian, the parameters refer to the constant term, amplitude 
    % term, stdev term and center term respectively.
    % For 2d circular gaussian: Parameters: p(1): z-offset, p(2): amplitude, 
    % p(3): xStdev, p(4): yStdev, p(5): xCenter, p(6): yCenter, p(7): tilt.
    
    if isempty(px) == 1 || isempty(py) == 1
        popt=[];
        resnorm=[];residual=[];ret=[];
        return
    end

    % First guess, without considering the tilt
    p0 = [((px(1)+py(1))/2) sqrt(abs(px(2)*py(2))) px(3) py(3) px(4) py(4)];
    
    % If tilt is specified, simply initialize the tilt with zero
    if tilt==1
        p0 = [p0 0];
    end

    %tic;
    options = optimset('Display','off',...
        'Algorithm','levenberg-marquardt','Jacobian','on','TolFun',1e-10,'TolX',1e-10);
    [popt,resnorm,residual,~,output] = lsqcurvefit(@gauss2dfunct,p0,grid,z,[],[],options);
    %time=toc;
    
    ret=output.iterations;
    
    %totalTime=toc;
end
