clearvars -except bead_nf bead_raw refZ refZ_fft

% Parameters
D0 = 30;
Di = 0;
thres = 1000;
nbors = 10;

%%

pxl_size = 54.2;
xax = 0:0.025:3;
xsurr_long = min(xax):0.01:max(xax);
parab = [];
sigma_xy=[];
totphoton=[];
% 
% pathname = 'C:\Users\cryst\OneDrive - University of Illinois - Urbana\Analysis\SELFI\Data\Test_221001\';
% pathname = 'C:\Users\cryst\OneDrive - University of Illinois - Urbana\Analysis\SELFI\Data\220614\';
pathname = 'C:\Users\cryst\OneDrive - University of Illinois - Urbana\Analysis\SELFI\Data\Test_220922\';
filename = 'test_AF647.tif';
FFT_tiff = [pathname 'MovieFFT_D0=' num2str(D0) '_' filename];
% filename = 'AF647-1nM_90pct641_21pct405_EMG100_200ms-crop.tif';
% filename = 'beads_4X_SELFI_0um_rep2.tif';
filename = 'Neuron_10nM_95pct641_10pct405_EMG50_100ms_HILO_SELFI_v5.tif';
% keyboard
[rawframes,numframe] = loadTiffMovie(pathname, filename);
%%
cf_all=rawframes;
% output_image=cell(1,numframe);
% %frames and peaks
% frpk = [];
% all_mols = {};
for NN = 4 : numframe
% for NN = 250:250
cf = cf_all{NN};
% FFT + LP filter
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
% keyboard
output_image{NN} = real((ifft2(double(G))));

% Peak finder
pks = FastPeakFind(output_image{NN},thres);
pks_pair = reshape(pks,2,[]); % reshapes pks into x-y coordinate pairs
dist = sqrt((pks(3:2:end-1)-pks(1:2:end-3)).^2 + (pks(4:2:end)-pks(2:2:end-2)).^2);
too_close = [find(dist<20) find(dist<20)+1];
pks_final = pks_pair; pks_final(:,too_close) = []; % remove pks too close together
[m,n] = size(pks_final);
curr_frpk = ones(n,3);
curr_frpk(:,1) = ones(n,1)*NN;
curr_frpk(:,2:3)=pks_final';
frpk = [frpk; curr_frpk];
% keyboard
% get individual interferograms
curr_fr_m = cell(n, 3); % 1=raw interferogram; 2=FT_img; 3=iFFT(G);
I0 = 5; I1 = 0;
sigg = [];
PF=[];
NP=[];
for ii = 1 : n
    if curr_frpk(ii,2) > 10 && curr_frpk(ii,2) < N-10 && curr_frpk(ii,3) > 10 && curr_frpk(ii,3) < M-10
        % Raw interferogram
        curr_fr_m{ii, 1} = cf(curr_frpk(ii,3)-nbors:curr_frpk(ii,3)+nbors, ...
            curr_frpk(ii,2)-nbors:curr_frpk(ii,2)+nbors);
        
        % FFT of raw interferogram, pre-lowpass
        [Mi , Ni] = size(curr_fr_m{ii ,1});
        FT_imgi = fft2(double(curr_fr_m{ii ,1}));
        curr_fr_m{ii, 2} = FT_imgi;
            
        % low-pass, then iFFT
        ui = 0:(Mi-1);
        idxi = find(ui>Mi/2);
        ui(idxi) = ui(idxi)-Mi;
        vi = 0:(Ni-1);
        idyi = find(vi>Ni/2);
        vi(idyi) = vi(idyi)-Ni;

        [Vi, Ui] = meshgrid(vi, ui);
        Dm = sqrt(Ui.^2 + Vi.^2);

        Hm = double( Dm <= I0 & Dm >=I1);
        Gm = Hm.*FT_imgi; 

        curr_fr_m{ii, 3} = real((ifft2(double(Gm))));
        
        %
        %Z axis localization w/ Pearson Coeff
        for zz= 1 : bead_nf
            R = corrcoef(imag(fftshift(Gm)),imag(fftshift(refZ_fft{zz})));
            PC(:,zz) = abs(R(2,1));
        end
%         xax_i=xax(68:92);
    PF(ii,:) = polyfit(xax,PC,2);
    [sig_x, sig_y, numphoton] = PSFfit(curr_fr_m{ii, 3},pxl_size);
    sigg(ii,:) = [sig_x,sig_y];
    NP(ii,:) = numphoton;
%     keyboard
    end
    
end
parab = [parab; PF];
sigma_xy = [sigma_xy;sigg];
totphoton = [totphoton;NP];
all_mols = vertcat(all_mols,curr_fr_m);
% keyboard
end