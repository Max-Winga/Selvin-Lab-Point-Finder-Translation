D0 = 3;
Di = 0;

pxl_size = 54.2;
xax = 1.6750:0.025:2.2750;
xsurr_long = min(xax):0.01:max(xax);

cf_all = frXY(~cellfun('isempty',frXY));
% cf_all = refZ;
sigma_xy = zeros(numel(cf_all),2 );
% 
pathname = 'C:\Users\cryst\OneDrive - University of Illinois - Urbana\Analysis\SELFI\Data\Test_221001\';
filename = 'test_AF647.tif';
FFT_tiff = [pathname 'FFT_D0=' num2str(D0) '_' filename];
% % N=29;
% 
% f = waitbar(0,'Please wait...');
%%
% hold on
for NN = 1 : numel(cf_all)
% for NN = 984:984
cf = cf_all{NN};

[M , N] = size(cf);
FT_img = fft2(double(cf));

u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

[V, U] = meshgrid(v, u);
D = sqrt(U.^2 + V.^2);

H = double( D <= D0 & D >=Di);
G = H.*FT_img;
% G = H.*test;

output_image(: , : , : , NN) = real(fftshift(ifft2(double(G))));
% output_image(:,:,:,NN) = abs(fftshift(double(G)));

% keyboard
% XY
% [sig_x, sig_y] = PSFfit(output_image(: , : , : , NN),pxl_size);
% sigma_xy(NN,1) = sig_x;
% sigma_xy(NN,2) = sig_y;
% keyboard
    for ii= 1 : bead_nf
        R = corrcoef(imag(fftshift(G)),imag(fftshift(refZ_fft{ii})));
        PC(:,ii) = abs(R(2,1));
        
    end

    Pearson_coef(NN,:) = PC(1:end);
%     PF(NN,:) = polyfit(xax(1:end),Pearson_coef(NN,:),2);
%     Pmax = find(polyval(PF(NN,:),xax) == max(polyval(PF(NN,:),xax)));
%     Pmax = find(Pearson_coef(NN,:)==max(Pearson_coef(NN,:)));
        PF(NN,:) = polyfit(xax,PC,2);
        
        if PF(NN,2)>0
            PF(NN,:)=PF(NN,:);
        elseif PF(NN,2)<0
            PF(NN,:)=-PF(NN,:);
        end
        Pmax_long = find(polyval(PF(NN,:),xsurr_long) == max(polyval(PF(NN,:),xsurr_long)));
%         Pmax_long = find(polyval(PF(NN,:),xsurr) == max(polyval(PF(NN,:),xsurr)));
%           keyboard  
     Zmax(NN) = xsurr_long(Pmax_long);
%     scatter(XYpos(NN,1),XYpos(NN,2),50, colorMap(Pmax_long,:),'filled');
% 
    if NN==1
        imwrite(uint16(output_image(:,:,:,NN)),FFT_tiff);
    else
        imwrite(uint16(output_image(:,:,:,NN)),FFT_tiff,'WriteMode','append');
    end

%     waitbar(NN/100)


end

% close(f);

%%
% subplot(1, 3, 1), imshow(cf, [min(min(cf)) max(max(cf))]), colorbar
% subplot(1, 3, 2), imshow(output_image(:,:,:,NN), [min(min(cf)) max(max(cf))]); colorbar
% subplot(1, 3, 3), pcolor(fftshift(H)); axis square tight; colorbar 