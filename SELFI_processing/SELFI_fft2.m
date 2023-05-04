
function [rf,ff,thetaf] = SELFI_fft2(rawframes)

    numframe = length(rawframes);
    
    for ii = 1 : numframe
        rf{ii} = real(fftshift(fft2(rawframes{ii})));
        ff{ii} = imag(fftshift(fft2(rawframes{ii})));
        thetaf{ii}= atand(ff{ii}./rf{ii});
        
    end

end