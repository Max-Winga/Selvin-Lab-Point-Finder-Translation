function [y,jac] = gauss1dfunct(p, xgrid)

  % For this function, we need to input the four parameters in p [p(1),
  % p(2), p(3), and p(4)], and the x range in data. The parameters refer to
  % the constant term, amplitude term, stdev term and center term
  % respectively.
  % The output will be y and the jacobian
  
  U0=(xgrid-p(4))/p(3);
  U=U0.*U0;
  U2=p(2)*U0/p(3);
  U3=U2.*U0;
  expU=exp(-U/2);
  
  y=p(1)+p(2)*expU;
  
  % The Jacobian of Gauss1d
  
  if nargout>1      % Jacobian required
      jac=zeros(length(xgrid),4);
      jac(:, 1)=1.0;
      jac(:, 2)=expU;
      jac(:, 3)=U3.*expU;
      jac(:, 4)=U2.*expU;
  end
  
end