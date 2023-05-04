% Run fitMolecule.m first

REF = INTERF_REF;

locs = zeros(numel(INTERF_final),4);

for ii = 1 : numel(INTERF_final)
    
    [sig_x, sig_y, zloc, polyf, N] = PSFfitWithZ(INTERF_final{ii},REF, pxl_size);
    locs(ii,:) = [sig_x, sig_y, zloc, N]; 
    fits{ii} = polyf;
end

