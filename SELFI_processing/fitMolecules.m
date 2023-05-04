pathname = 'D:\Data\SELFI_230215\';
addpath(pathname);
filename = 'TSB_488_withPH_p400nm.tif';
% filename = 'svaf2nM_50ms_100EMG_640nm0p6_m200nm.tif';
% keyboard
[rawframes,numframe] = loadTiffMovie(pathname, filename);
curr_frame = rawframes{1};
%%
pxl_size = 43.3;

nbors = 15;
thres = 1250;
% thres = 1.5e4;
D0 = 25;
D1 = 1;
R = 20;

Z_uncertainty = [];
XYpos = [];
xy_uncertainty = [];
cent_final = zeros(2,1);
MOLE_final = cell(1,1);
INTERF_final = cell(1,1);
for ii = 1 :  numframe
    curr_frame = rawframes{ii};
    [B1,B2]=Gaussian_image_filtering(curr_frame,D0);
    gg = B1;
    [L,W] = size(gg);
    subplot(1,2,1); 
    [cent,varargout]=FastPeakFind(gg,thres);
    cent=reshape(cent,2,[]);
    cent(:,cent(1,:)<=nbors | cent(1,:)>=W-nbors)=[];
    cent(:,cent(2,:)<=nbors | cent(2,:)>=L-nbors)=[];
    [~,numLoc]= size(cent);

    [xx,yy] = meshgrid(1:L,1:W);

    % BUILD DISTANCE MATRIX TO REMOVE MOLECULES THAT ARE TOO CLOSE TOGETHER
    for jj = 1:numLoc
%         dist = sqrt(( cent(1,jj+1)-cent(1,jj)).^2  + (cent(2,jj+1)-cent(2,jj)).^2);
%         if floor(dist) > nbors
            
%             MOLE = double(gg(cent(2,jj)-nbors:cent(2,jj)+nbors, cent(1,jj)-nbors:cent(1,jj)+nbors));
%             INTERF = double(curr_frame(cent(2,jj)-nbors:cent(2,jj)+nbors, cent(1,jj)-nbors:cent(1,jj)+nbors));
            R2C = (xx - cent(1,jj)).^2/R^2 + (yy-cent(2,jj)).^2/R^2;
            circlefilt = (R2C < 1); %circular mask
            
            ISOLATE= gg.*circlefilt; 
%             MOLE = double(ISOLATE(cent(2,jj)-nbors:cent(2,jj)+nbors, cent(1,jj)-nbors:cent(1,jj)+nbors));
%             MOLE(MOLE==0) = MOLE(15,2);
            cISO = double(curr_frame).*circlefilt;
            INTERF = double(cISO(cent(2,jj)-nbors:cent(2,jj)+nbors, cent(1,jj)-nbors:cent(1,jj)+nbors));
%             INTERF(INTERF==0) = INTERF(15,2)
% keyboard
            [MOLE,~]=Gaussian_image_filtering(INTERF,D1);
%             keyboard
            [sig_x, sig_y, N] = PSFfit(abs(MOLE), pxl_size);
            xy_uncertainty = [xy_uncertainty;[sig_x, sig_y, N]];

            if isempty(N) == 0
                cent_final = [cent_final,cent(:,jj)];
                MOLE_final{end+1} = MOLE; 
                
                INTERF_final{end+1} = INTERF;
            end
            
%         end
    end
%     subplot(1,2,2);
    cent_final = cent_final(:,2:end);
%     plotIm(curr_frame); hold on; scatter(cent_final(1,:),cent_final(2,:),80,'ro');caxis([800 3e3])
end
MOLE_final=MOLE_final(~cellfun('isempty',MOLE_final));
INTERF_final=INTERF_final(~cellfun('isempty',INTERF_final));

distFromCenter = sqrt((cent_final(1,:)-512).^2+(cent_final(2,:)-512).^2);
mindist = find(distFromCenter==min(distFromCenter));

% t1=intersect(find(distFromCenter<100),find(xy_uncertainty(:,3)<2e6));
% t2=intersect(find(distFromCenter<200),find(xy_uncertainty(:,3)<2e6));
% t3=intersect(find(distFromCenter<300),find(xy_uncertainty(:,3)<2e6));