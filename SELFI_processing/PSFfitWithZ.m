% IMAGE=output_image(:,:,:,NN);
% IMAGE=double(cf);

function [sig_x, sig_y, zloc, polyf, N] = PSFfitWithZ(IMAGE,REF, pxl_size)


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

a = pxl_size;

bkg = reshape(IMAGE-z,1,[]);
b = std(bkg);

sig_ax2 = si_x^2 + a^2/12;
sig_ay2 = si_y^2 + a^2/12;

sig_x = sqrt(sig_ax2/N * (16/9 + 8*pi*sig_ax2*b^2/(N*a^2)));
sig_y = sqrt(sig_ay2/N * (16/9 + 8*pi*sig_ay2*b^2/(N*a^2)));


% Z localization
[PSpec,~] = powerSpec(IMAGE);
PCoeffCorr = zeros(1, numel(REF));

for ii = 1 : numel(REF)
    [REF_PSpec,~] = powerSpec(REF{ii});
    RR = corrcoef(REF_PSpec, PSpec);
    PCoeffCorr(ii) = RR(1,2);
end

REF_ax = -floor(numel(REF)/2):1:floor(numel(REF)/2);
polyf = fit(REF_ax', PCoeffCorr','poly2');
zloc = roots(polyder(coeffvalues(polyf)));
% keyboard

else
    sig_x=0; sig_y=0; zloc = 0; N=0; polyf = 0;
end
%%
end