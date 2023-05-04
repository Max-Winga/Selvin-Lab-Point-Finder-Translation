

function [rawframes,numframe] = loadTiffMovie(pathname, filename)
% keyboard
ftbo = [pathname filename];
info = imfinfo(ftbo);
numframe = length(info);
rawframes={};
for ii = 1 : numframe
    rawframes{ii} = imread(ftbo, ii);
    
end

% cookedframes = mat2gray(rawframes);
% implay(cookedframes);
end