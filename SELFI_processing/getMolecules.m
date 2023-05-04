pathname = 'C:\Users\wjchen10\OneDrive - University of Illinois - Urbana\Analysis\SELFI\Data\SELFI_221024\';
addpath(pathname);
% filename = 'AF647-1nM_95pct641_10pct405_EMG50_100ms_10269_newbuffer_HILO.tif';
filename = 'Sil60X_X10269_640P100_1nMsvAF_1-crop.tif';
% keyboard
[rawframes,numframe] = loadTiffMovie(pathname, filename);

% X1Y1 = 'AF647-1nM_95pct641_10pct405_EMG50_100ms_10269_newbuffer_HILO-allfr_XYZ.csv';
% X1Y1 = 'AF647-1nM_90pct641_21pct405_EMG100_200ms-allfr_XYZ.csv';
% opts = detectImportOptions(X1Y1);
% opts.SelectedVariableNames = [2:5,9];
% frXY_dex = readmatrix(X1Y1, opts);
% frXY={};

pxl_size = 216.7;

nbors = 10 ;
Z_uncertainty = [];
XYpos = [];
% isolate frames
%% 

for ii = 1 : max(frXY_dex(:,1))
    curr_fr = frXY_dex(frXY_dex(:,1) == ii, :);
    if isempty(curr_fr)==0
    curr_fr_m = cell(1, numel(curr_fr(:,1)));
    curr_raw = rawframes(:,:,:,ii);
        for jj = 1 : numel(curr_fr(:,1))
            Xc = round(curr_fr(jj, 3)/pxl_size);
            Yc = round(curr_fr(jj, 2)/pxl_size); %keyboard
            if Xc > 10 && Xc < 1014 && Yc > 10 && Yc < 1014
                curr_fr_m{jj} = curr_raw(Xc-nbors:Xc+nbors, Yc-nbors:Yc+nbors);
                frXY{ii,jj} = curr_fr_m{jj};
                XYpos = [XYpos; Xc Yc];
%                 Pearson_coef = zeros(numel(curr_fr_m),znum-20);
%                 
%                 for kk= 1 : znum-21
%                     R = corrcoef(double(curr_fr_m{jj}),double(refZ(:,:,:,kk)));
%                     PC(:,kk) = R(2,1);
%                     
%                 end
%                 
%                 Pearson_coef(jj,:) = PC(21:end);
%                 PF(jj,:) = polyfit(xax,Pearson_coef(jj,:),2);
%                 Pmax = find(polyval(PF(jj,:),xax) == max(polyval(PF(jj,:),xax)));
%                 Zmax(jj) = xax(Pmax);
                
            end
        end
%     Z_uncertainty = [Z_uncertainty, Zmax];
    else
        continue
    end
end