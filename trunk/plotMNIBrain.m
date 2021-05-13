function [H] = plotMNIBrain(alpha,bg)
% PLOTMNIBRAIN plots the freesurfer MNI surface on a black background. 
%   
%   plotMNIbrain(alpha) uses the amount of transparency in alpha for the 
%   MNI surface. 
% 
%   plotMNIbrain(alpha,bg) allows the user to specify the color of the
%   background as 'w' (white) or 'k' (black).
% 
%   wrapper for the plotFreeSurf function, which is a wrapper for some
%   freesurfer functions. 
% 

% author EHS20160317 

surfPath = '/home/user1/code/matlab/brainVisualization/trunk/MNI_fsl_surf/';
surfName = 'pial';

if ~exist('bg','var') || strcmp(bg,'k')
    plotFreeSurf(surfPath,surfName,'both',alpha)
elseif strcmp(bg,'w')
    plotFreeSurf_wbg(surfPath,surfName,'both',alpha)
else
    [V,F,p,H] = plotFreeSurf_bbg(surfPath,surfName,'both',alpha)
end

end
