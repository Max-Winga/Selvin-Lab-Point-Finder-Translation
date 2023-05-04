beadpath = 'C:\Users\wjchen10\OneDrive - University of Illinois - Urbana\Desktop\LabView_Test\Test_221011\';
beadfile = 'Beads_SCMOS_488nm30pct_adjust-1st-lens-X773Y1253.tif';
PSF_movie = 'bead_PSF_ifft2.tif';
FFT_tiff = [beadpath PSF_movie];

[bead_raw,bead_nf] = loadTiffMovie(beadpath, beadfile);
%%
% nbors = 10;
% XZs = 81; 
% YZs = 592;
% refZ = {};
% refZ_fft={};
D0=3;
Di=0;

testz_targ={};
    
% isolate frames
for ii = 1 :  bead_nf
    curr_bead = bead_raw{ii};
% keyboard
    refZ{ii} = uint16(curr_bead);
    [M , N] = size(refZ{ii});
    FT_img = fft2(double(refZ{ii}));
    
    u = 0:(M-1);
    idx = find(u>M/2);
    u(idx) = u(idx)-M;
    v = 0:(N-1);
    idy = find(v>N/2);
    v(idy) = v(idy)-N;
    
    [V, U] = meshgrid(v, u);
    D = sqrt(U.^2 + V.^2);
    
    H = double( D <= D0 & D >=Di);
    refZ_fft{ii} = FT_img;
    testz_targ{ii,1} = double(refZ{ii});
    testz_targ{ii,2} = imag(fftshift(refZ_fft{ii}));
    
    % G = H.*test;
%     keyboard
    % output_image(: , : , : , NN) = real(ifft2(double(G)));
    towrite = ifft2(refZ_fft{ii});
    testz_targ{ii,3}=towrite;
    if ii==1
        imwrite(uint16(towrite),FFT_tiff);
    else
        imwrite(uint16(towrite),FFT_tiff,'WriteMode','append');
    end
end

% clearvars -except bead_nf bead_raw refZ refZ_fft testz