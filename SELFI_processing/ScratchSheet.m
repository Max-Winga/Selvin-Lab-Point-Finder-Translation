%%
fftx = zeros(2^nextpow2(M+2), 2^nextpow2(N+2));
ffty = zeros(2^nextpow2(M+2), 2^nextpow2(N+2));
for ii=1:176
    fftx(ii,:)=(fft(p200(ii,:)-mean(p200(ii,:)),2^nextpow2(M+2)));
    ffty(:,ii)=(fft(p200(:,ii)-mean(p200(:,ii)),2^nextpow2(M+2)));
end


%% selfi theory testing
lambda=  550/1e3; %um
NA = 1.4;
r0 = 1.22*lambda/(2*NA);
k = 2*pi/lambda; %um^-1
I0 = 1;

xy = 0.25;
deltaxy = 0.001;
rr = linspace(-xy,xy,2^9);
[X,Y] = meshgrid(rr,rr);
% r = sqrt(X.^2+Y.^2);
[theta,r] = cart2pol(X,Y);
z = (-1:0.1:1); %um
sig0 = r0./(2*sqrt(2*log(2)));

R = z.*(1+4*k^2*sig0^4./z.^2);

r = r./lambda;

[L,W]=size(r);
phi = zeros(L,W,length(R));
PSF = zeros(L,W,length(R));
SIG = zeros(1,length(R));
for ii = 1:length(R)
    SIG(ii) = sig0*(1+z(ii).^2./(8*k^2*sig0^4));
    phi(:,:,ii) = -k*(r.^2-r0.^2)./R(ii);
    PSF(:,:,ii) = I0./(2*pi*SIG(ii).^2).*exp(-r.^2./(2*SIG(ii).^2)).*exp(-1i*phi(:,:,ii));
end

%% selfi testing
rectang_func = rectang(2^8,rr);
IFT = fftshift(fft2(PSF(:,:,10)));
interferogram = IFT.*rectang_func;
imagesc(abs(interferogram));
axis equal tight
zoom on