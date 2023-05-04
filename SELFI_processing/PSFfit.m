% IMAGE=output_image(:,:,:,NN);
% IMAGE=double(cf);

function [sig_x, sig_y, N] = PSFfit(IMAGE,pxl_size)


[ny,nx]=size(abs(IMAGE));
grid=[nx,ny,1:nx,1:ny];
% keyboard
[popt,~,~,~]=gauss2dfit(abs(IMAGE), grid,0);
% keyboard
if isempty(popt) == 0

xp=reshape(repmat((grid(3:nx+2)-popt(5))/popt(3),[ny 1]),ny*nx,1);       % Expand x values, make 1d
yp=reshape(repmat((grid(nx+3:end)'-popt(6))/popt(4),[1 nx]),ny*nx,1);    % Expand y values, make 1d
expU=exp(-0.5*(xp.*xp+yp.*yp));     % Find exp term
z=popt(1)+popt(2)*expU; 

N = sum(z-popt(1));
% keyboard
z=reshape(z,[ny,nx]);

si_x = popt(3)*pxl_size;
si_y = popt(4)*pxl_size;
% si_x = 141.02;
% si_y = 141.02;
a = pxl_size;

bkg = reshape(IMAGE-z,1,[]);
b = std(bkg);

% sig_x = sqrt(si_x.^2./N + a^2/12./N + 8*pi*si_x^4*b^2/a^2./N.^2);
% sig_y = sqrt(si_y.^2./N + a^2/12./N + 8*pi*si_y^4*b^2/a^2./N.^2);

sig_ax2 = si_x^2 + a^2/12;
sig_ay2 = si_y^2 + a^2/12;

sig_x = sqrt(sig_ax2/N * (16/9 + 8*pi*sig_ax2*b^2/(N*a^2)));
sig_y = sqrt(sig_ay2/N * (16/9 + 8*pi*sig_ay2*b^2/(N*a^2)));
% keyboard
else
    sig_x=[]; sig_y=[]; N=[];
end
%%
end