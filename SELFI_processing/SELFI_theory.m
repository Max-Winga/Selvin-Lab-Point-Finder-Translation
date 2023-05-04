grating = zeros(L,W);
p = 20; %phase mask period, um
for nn = -50:50
    for mm= -50:50
        grating = grating + sinc(pi/2*(2*nn+1))*sinc(pi/2*(2*mm+1))*exp(-2*1i*pi/p*((2*nn+1)*X + (2*mm+1)*Y)); 
%         pcolor(abs(grating)); shading flat
    end
end