function [z,jac]=gauss2dfunct(p,grid)

% 2d gaussian function that allows calculation for elliptical gaussian,
% with certain tilt. This gaussian function is adapted from the IDL version
% of GAUSS2_FUNCT.

% Input: 
% 'p' is an array of 6 or 7 elements. If it has 7 elements, tilt will
% be included. The parameters are: z-offset, amplitude, x-stdev, y-stdev,
% x-center, y-center, tilt.
% 'grid' is an array consisting of the dimension in x (grid(1)) and y (grid(2)),
% and the gridding in x and y (which will be placed right after the 
% dimensions, in the same column array). This input is [nx ny x y].

% Output:
% The outputs are the 2d gaussian function denoted as 'z', and the jacobian
% denoted as 'jac'.

nx=grid(1);             % Obtain the dimensions in x and y
ny=grid(2);

% Jacobian of 2d gaussian function

m=length(p);            % Find out the number of parameters
jac=zeros(ny*nx,m);     % Allocate memory for jac

tilt=length(p)==7;
if tilt==1;
    xpt=repmat((grid(3:nx+2)-p(5)),[ny 1]);         % Expand x values
    ypt=repmat((grid(nx+3:end)'-p(6)),[1 nx]);      % Expand y values
    s=sin(p(7));                                    % Values of sin and cos
    c=cos(p(7));
    xp=reshape(xpt*(c/p(3)) - ypt*(s/p(3)),ny*nx,1);   % Make 1d
    yp=reshape(xpt*(s/p(4)) + ypt*(c/p(4)),ny*nx,1);   % Make 1d
    expU=exp(-0.5*(xp.*xp+yp.*yp));     % Find exp term
    z=p(1)+p(2)*expU;                   % Calculate z
    
    % Determine partial derivatives
    if nargout>1        % Jacobian required
        pexpU=p(2)*expU;
        Ux=pexpU.*xp/p(3);
        Uy=pexpU.*yp/p(4);
        jac(:,1)=1;
        jac(:,2)=expU;
        jac(:,3)=Ux/p(3);
        jac(:,4)=Uy/p(4);
        jac(:,5)=Ux*c+Uy*s;
        jac(:,6)=-Ux*s+Uy*c;
        jac(:,7)=-pexpU.*xp.*yp*(p(3)/p(4)-p(4)/p(3));
    end
    
else
    xp=reshape(repmat((grid(3:nx+2)-p(5))/p(3),[ny 1]),ny*nx,1);       % Expand x values, make 1d
    yp=reshape(repmat((grid(nx+3:end)'-p(6))/p(4),[1 nx]),ny*nx,1);    % Expand y values, make 1d
    expU=exp(-0.5*(xp.*xp+yp.*yp));     % Find exp term
    z=p(1)+p(2)*expU;                   % Calculate z
    
    % Determine partial derivatives
    if nargout>1        % Jacobian required
        pexpU=p(2)*expU;
        Ux=pexpU.*xp/p(3);
        Uy=pexpU.*yp/p(4);
        jac(:,1)=1;
        jac(:,2)=expU;
        jac(:,3)=Ux/p(3);
        jac(:,4)=Uy/p(4);
        jac(:,5)=Ux;
        jac(:,6)=Uy;
    end
end

end
