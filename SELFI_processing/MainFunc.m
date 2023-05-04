pathname = 'C:\Users\wjchen10\OneDrive - University of Illinois - Urbana\Analysis\SELFI\Data\Test_220922\';
filename = 'Neuron_10nM_95pct641_10pct405_EMG50_100ms_HILO_SELFI_v4.tif';
Zname = 'TSB_10pct641_EMG50_100ms_zstep25nm_0-3um.tif';

[ZStack, znum] = loadTiffMovie(pathname, Zname);

%%
Pearson_coef = zeros(numel(curr_fr_m),znum-20);

for jj = 1 : numel(curr_fr_m)
    for ii= 1 : znum-21
        R = corrcoef(double(curr_fr_m{jj}),double(refZ(:,:,:,ii)));
        PC(:,ii) = R(2,1);
        
    end
    
    Pearson_coef(jj,:) = PC(21:end);
    PF(jj,:) = polyfit(xax,Pearson_coef(jj,:),2);
    Pmax = find(polyval(PF(jj,:),xax) == max(polyval(PF(jj,:),xax)));
    Zmax(jj) = xax(Pmax);
end
