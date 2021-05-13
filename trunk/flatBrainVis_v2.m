% For working with Tyler's much faster (as usual) GUI.
surfPath = '/home/user1/Desktop/R5/Tyler/CorticalStim/CS201810/Imaging/Registered/';
surfTrodesName = 'electrodes.mat';
alpha = 0.7;
patientID = 'u201810';

% loading data
load(fullfile(surfPath,surfTrodesName))

% changing variable names to match the rest of the script.
% [20181030] TODO:: fix this to use the vraiable names from Tyler's GUI.
CentroidCoordinates = ElecXYZProjRaw;
ncoords = length(CentroidCoordinates);
F = BrainSurfRaw.faces;
V = BrainSurfRaw.vertices;

% if you just want to visualize electrodes on the pial surface.
plotOGTrodes = false;
if plotOGTrodes
    figure('Color',rgb('black'));
    % plotting the brain surface
    hold ongcc
    p=patch('faces',F,'vertices',V, 'facecolor', 'flat',  'edgecolor', 'none', 'facealpha', alpha);
    scatter3(CentroidCoordinates(:,1),CentroidCoordinates(:,2),CentroidCoordinates(:,3),30,rgb('lime'),'filled');
    hold off
    
    % set face color
    faceColor = repmat([1 1 1], length(F), 1);
    set(p,'FaceVertexCData',faceColor);
    
    % adjust view.
    daspect([1 1 1])
    view(3); axis tight off
    view([50 -40 100])
    camlight
    lighting gouraud
end


%% finding nearest vertices within a threshold for each hemisphere.
hsFlag = 'both';  % can be 'left', 'right', 'both'
switch hsFlag
    % is it still this simple?
    case {'left'}
        hemisphereName = 'leftHemi';
        ContactIdcs = find(CentroidCoordinates(:,1)<0);
        vertexIdcs = find(V(:,1)<0);
    case {'right'}
        hemisphereName = 'rightHemi';
        ContactIdcs = find(CentroidCoordinates(:,1)>0);
        vertexIdcs = find(V(:,1)>0);
    case {'both'}
        hemisphereName = 'bothHemis';
        ContactIdcs = 1:length(CentroidCoordinates(:,1));
        vertexIdcs = 1:length(V(:,1));
end
% faces on a particular hemisphere.
[FVintersect,faceIdcs,~] = intersect(F(:,1),vertexIdcs);

% if you just want to visualize electrodes on the pial surface
% (separated by hemisphere)
plotHemiTrodes = true;
if plotHemiTrodes
    figure('Color',rgb('black'));
    % plotting the brain surface
    hold on
    
    %         p=patch('faces',F,'vertices',V, 'facecolor', 'flat',  'edgecolor', 'none', 'facealpha', alpha);
    %                 % set face color
    %         faceColor = repmat([1 1 1], length(F), 1);
    %         set(p,'FaceVertexCData',faceColor);
    
    scatter3(V(vertexIdcs,1),V(vertexIdcs,2),V(vertexIdcs,3),2,rgb('silver'),'filled');
    scatter3(CentroidCoordinates(ContactIdcs,1),CentroidCoordinates(ContactIdcs,2),CentroidCoordinates(ContactIdcs,3),30,rgb('crimson'),'filled');
    hold off
    
    
    % axis details
    daspect([1 1 1])
    view(3); axis tight off
    camlight
    lighting gouraud
end

% finding the closest vertex to each electrode.
for el = 1:length(ContactIdcs)
    [minDist(el),Vinds(el)] = min(sqrt((CentroidCoordinates(ContactIdcs(el),1)-V(:,1)).^2+(CentroidCoordinates(ContactIdcs(el),2)-V(:,2)).^2+(CentroidCoordinates(ContactIdcs(el),3)-V(:,3)).^2));
end

% [20181030]: I think the coordinates are now in millimeters.
%~~~~THRESHOLD~~~~
distThresh = 7;
%~~~~~~~~~~~~~~~~~

% discarding distant electrodes.
elInds = minDist<distThresh;
sEEGvertices = V(Vinds(elInds),:);
adjacentFaces = {};

% if you want to project the electrodes on the flat surface
% mostly for configuration.
projectTrodes = true;
if projectTrodes
    hold on
    pts = scatter3(sEEGvertices(:,1),sEEGvertices(:,2),sEEGvertices(:,3),20,rgb('goldenrod'),'filled')
    hold off
end

plotHistogram = true;
if plotHistogram
    figure
    histogram(minDist,'facecolor','k')
    xlabel('mri distance (mm)')
    ylabel('count')
end

% ^^^^^^^Electrode Stuff^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% [20181030] surface dimensionality reduction (attempts) below.

% NOTES: I guess I want to take the vertices that are closest to the
% 3-D plane Y = Z and delete those, then make automated cuts on the
% medial surface of the brain, then use isomap to dimensionality
% reduce the vertices from three dimensions to two dimensions. I'm
% thinkin gof discarding the face information at this point for
% simplicity's sake. I'm not sure how to deal with what happens to a
% face when a single vertex is removed.

% [20181105] going to just try isomap without the cuts, since the
% cerebellum is already there. Doing the Isomapping from
% www.numerical-tours.com/matlab/meshdeform_3_flattening
% -- this failed.

% plotting 2-D surface.
twoDhemi = false;
if twoDhemi
    % [20181105] now trying locally linear embedding
    nNeibs = 6;
    Y = lle(V(vertexIdcs,:)',nNeibs,2);
    
    figure('Color',rgb('black'));
    hold on
    scatter(Y(1,:),Y(2,:),3,rgb('silver'),'filled')
    
    hold off
    axis off
end
% this looks crazy -- definitely a failure.


%% loading seizure and calculating phase-locked high gamma
seizurePath = '/home/user1/data/Seizures/UofU';
seizureFileName = '201810.edf';
dataFormat = 'edf';
trimLims = [1e5 5e5];
if isequal(dataFormat,'blackrock')
%     openNSx('/media/user1/data4TB/Seizures/transferredData/CUBF22/20170320-162443-038.ns3')
    labels = deblank({NS3.ElectrodesInfo.Label});
    labels = labels(1:92);
    tmp = double(NS3.Data(1:92,:));
    Fs = NS3.MetaTags.SamplingFreq;
    
    dataRange = [4e5 7.5e5]; % range of data to save in samples.
    for ch = 1:length(labels)
        data(ch,:) = resample(tmp(ch,dataRange(1):dataRange(2)),500,Fs);
    end
    Fs = 500;
elseif isequal(dataFormat,'edf')
    % loading data
    [hdr,data] = edfread(fullfile(seizurePath,seizureFileName));
    Fs = 500;
    % demeaning and trimming data over channels and time. 
    tmp = data(1:length(CentroidCoordinates),trimLims(1):trimLims(2))-mean(data(1:length(CentroidCoordinates),1:trimLims(1)),2);
    data = tmp;
    tSec = linspace(0,floor(size(data,2)/Fs),size(data,2));
else
    fprintf('\nStill working on this...')
end


%% dominant Frequency analyses.
dFname = sprintf('/media/user1/data4TB/Seizures/transferredData/%s/%s_dominantFrequencyData.mat',patientID,patientID);
if ~exist(dFname)
    [dF] = dominantFrequency(patientID,data,Fs,[],false);
    save(dFname,'dF','-v7.3')
else
    fprintf('you already calculated dominant frequency for %s seizure. loading those...',patientID)
    load(dFname,'dF')
end
nChans = length(CentroidCoordinates);


%% phase locked high gamma calculation over time
% for each discharge, caluclate the PLHG for each channel.
phiLo = nanmean(dF.PHIft(1:nChans,dF.fHz>4 & dF.fHz<30,:),2);
phiHi = nanmean(dF.PHIft(1:nCHans,dF.fHz>80 & dF.fHz<150,:),2);

% which signal to phase lock
hGamSig = nanmean(dF.Sft(1:nChans,dF.fHz>80 & dF.fHz<=150,:),2);
lGamSig = nanmean(dF.Sft(1:nChans,dF.fHz>30 & dF.fHz<=50,:),2);

% phase locked signal
PLHG = abs(squeeze(hGamSig.*(phiLo-phiHi)));
PLHGn = 100*round((PLHG-min(PLHG,[],2))./max(PLHG-min(PLHG,[],2),[],2),2);

% smoothing and normalizing
PLHGs = smoothdata(PLHG,2,'gaussian',2*Fs);
PLHGsn = 100*round((PLHGs-min(PLHGs,[],2))./max(PLHGs-min(PLHGs,[],2),[],2),2);

visualizeSeizureData = true;
if visualizeSeizureData
    figure
    subplot(2,2,1)
    imagesc(hGamSig)
    axis square tight
    subplot(2,2,2)
    imagesc(PLHGn)
    axis square tight
    subplot(2,2,3)
    imagesc(PLHGsn)
    axis square tight
end


%% plotting the PLHG across faces on the inflated brain
% heat map showing percentage of PLHG activation.
heatMap = flipud(colormap(hot(100)));
alpha = 0.3;

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% using vertices in visualization
% number of nearest faces to use.
kVerts = 100;
[VindExpansion,D] = knnsearch(V,[Vinds' Vinds' Vinds'],'K',kVerts);

% converting electrodes to faces again.
PLHGexpansion = zeros(0,size(PLHGsn,2));
for fc = 1:length(Vinds)
    % implements colormap smoothing
    catVar = uint8(repmat(PLHGsn(fc,:),kVerts,1).*repmat(linspace(1,0.25,kVerts)',1,size(PLHGsn,2)));
    % builds matrix of faces
    PLHGexpansion = cat(1,PLHGexpansion,catVar);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


keyboard
%% plotting the images and making the movie...
% movie params
introFrames = 11;
origAz = 270;
origEl = 0;
propAz = 220;
propEl = -10;
moveAz = linspace(origAz,propAz,introFrames);
moveEl = linspace(origEl,propEl,introFrames);
nFrames = size(PLHGfaces,2);

% now making the movie.
for frm = [1:introFrames introFrames+1:250:nFrames]
    fprintf('\nplotting frame %d of %d...',frm,nFrames)
    % now looping over frames
    if frm==1
        % plotting surface once.
        figure('Color',rgb('black'));
        % plotting the brain surface
        p=patch('faces',F,'vertices',V, 'facecolor', 'flat',  'edgecolor', 'none', 'facealpha', alpha);
        
        % set face color
        faceColor = repmat([1 1 1], length(F), 1);
        set(p,'FaceVertexCData',faceColor);
        
        % axis details
        daspect([1 1 1])
        axis tight off
        
        % setting up movie stuff
        dbPath = '/media/user1/data4TB/Dropbox/Dropbox/';
        movie_file_name = [dbPath patientID 'flattenedBrain_4xspeed.avi'];
        movie_object = VideoWriter(movie_file_name,'Motion JPEG AVI');
        movie_object.FrameRate = 10;
        movie_object.Quality = 50;
        open(movie_object);
        
        % setting up figure
        movie_handle = figure(gcf);
        
        % vis deets
        halfMaximize(gcf,'left')
        view(origAz,origEl)
        
        camlight
        lighting gouraud
        
        % write frame
        writeVideo(movie_object,getframe(movie_handle));
        
    elseif frm<introFrames %this implements a move in the view at first, so that
        alpha = 1/frm;
        % plotting a disappearing surface and the scatter
        hold on
        p=patch('faces',F,'vertices',V, 'facecolor', 'flat',  'edgecolor', 'none', 'facealpha', alpha);
        scatter3(V(vertexIdcs,1),V(vertexIdcs,2),V(vertexIdcs,3),2,rgb('silver'),'filled');
        scatter3(CentroidCoordinates(ContactIdcs,1),CentroidCoordinates(ContactIdcs,2),CentroidCoordinates(ContactIdcs,3),30,rgb('crimson'),'filled');
        hold off
        
        % slight chagne in camera view
        view(moveAz(frm-1),moveEl(frm-1))
        
        % write frame
        writeVideo(movie_object,getframe(movie_handle));
        
        % adjust view for next section.
        view(propAz,origEl)
        
    elseif frm>=introFrames
        %         try
        % selectively coloring faces.
        faceMap = zeros(size(PLHGfaces,1),3);
        for fm=1:size(PLHGfaces,1)
            faceMap(fm,:) = heatMap(PLHGfaces(fm,frm),:);
        end
        faceColor(VindExpansion,:) = faceMap;
        set(p,'FaceVertexCData',faceColor);
        
        % write frame
        writeVideo(movie_object,getframe(movie_handle));
        
        % removing plot
        faceColor = repmat([1 1 1], length(F), 1);
        set(p,'FaceVertexCData',faceColor);
        %         catch
        %             fprintf('\nunable to plot frame\n')
        %             close(movie_object)
        %         end
        
    end
end
close(movie_object)


