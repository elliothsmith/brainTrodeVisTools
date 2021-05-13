% params
surfPath = '/usr/local/freesurfer/subjects/CUBF22/surf/';
alpha = 1;
patientID = 'CUBF22';


%% loading centroids
load('/media/user1/data4TB/Dropbox/Dropbox/CUBF22_and_26/CUBF22_LOC/CUBF22FinalCoords_std.mat')
CentroidCoordinates = coords_new;
ncoords = length(CentroidCoordinates);
lhContacts = find(CentroidCoordinates(:,1)<0);
rhContacts = CentroidCoordinates(:,1)>0;


%% load convoluted pial Surface
% the vertices and faces correspond between pial and inflated surfaces.
[V,F] = freesurfer_read_surf([surfPath 'lh.pial']);

% if you just want to visualize electrodes on the pial surface.
plotOGTrodes = false;
if plotOGTrodes
    [V,F,~] = plotFreeSurf_bbg(surfPath,'pial','lh',0.4);
    hold on
    scatter3(CentroidCoordinates(lhContacts,1),CentroidCoordinates(lhContacts,2),CentroidCoordinates(lhContacts,3),30,rgb('springgreen'),'filled');
    hold off
    view(200,30)
    saveas(gcf,'./CUBF22_OGtrodes.pdf')
end


%% finding nearest vertices within a threshold
for el = 1:length(lhContacts)
    [minDist(el),Vinds(el)] = min(sqrt((CentroidCoordinates(lhContacts(el),1)-V(:,1)).^2+(CentroidCoordinates(lhContacts(el),2)-V(:,2)).^2+(CentroidCoordinates(lhContacts(el),3)-V(:,3)).^2));
end

plotHistogram = false;
if plotHistogram
    figure
    histogram(minDist, 50,'facecolor','k')
    xlabel('mri distance')
    ylabel('count')
end

%~~~~THRESHOLD~~~~
distThresh = 10;
%~~~~~~~~~~~~~~~~~

% discarding distant electrodes.
elInds = minDist<distThresh;
sEEGvertices = V(Vinds(elInds),:);
adjacentFaces = {};

% if you want to project the electrodes on the flat surface
% mostly for configuration.
projectTrodes = false;
if projectTrodes
    hold on
    pts = scatter3(sEEGvertices(:,1),sEEGvertices(:,2),sEEGvertices(:,3),20,[0 0.2 0.3],'filled')
    hold off
end

% either saving or closing the flat surface.
% TODO:: replace the
visFlatBrain=false;
if visFlatBrain
    % vis deets
    halfMaximize(gcf,'left')
    view(200,0)
    saveas(gcf,'/media/user1/data4TB/Dropbox/Dropbox/flatBrainContacts.pdf')
    close(gcf)
end


%% loading seizure and calculating phase-locked high gamma
openNSx('/media/user1/data4TB/Seizures/transferredData/CUBF22/20170320-162443-038.ns3')
labels = deblank({NS3.ElectrodesInfo.Label});
labels = labels(1:92);
tmp = double(NS3.Data(1:92,:));
Fs = NS3.MetaTags.SamplingFreq;

dataRange = [4e5 7.5e5]; % range of data to save in samples.
for ch = 1:length(labels)
    data(ch,:) = resample(tmp(ch,dataRange(1):dataRange(2)),500,Fs);
end
Fs = 500;


%% dominant Frequency analyses.
dFname = sprintf('/media/user1/data4TB/Seizures/transferredData/CUBF22/%s_dominantFrequencyData.mat',patientID);
if ~exist(dFname)
    [dF] = dominantFrequency(patientID,data,Fs,[],false);
    save(dFname,'dF','-v7.3')
else
    fprintf('you already calculated dominant frequency for %s seizure. loading those...',patientID)
    load(dFname,'dF')
end


%% phase locked high gamma calculation over time
% for each discharge, caluclate the PLHG for each channel.
phiLo = nanmean(dF.PHIft(:,dF.fHz>4 & dF.fHz<30,:),2);
phiHi = nanmean(dF.PHIft(:,dF.fHz>80 & dF.fHz<150,:),2);

% which signal to phase lock
hGamSig = nanmean(dF.Sft(:,dF.fHz>80 & dF.fHz<=150,:),2);
lGamSig = nanmean(dF.Sft(:,dF.fHz>30 & dF.fHz<=50,:),2);

% phase locked signal
PLHG = abs(squeeze(hGamSig.*(phiLo-phiHi)));
PLHGn = 100*round((PLHG-min(PLHG,[],2))./max(PLHG-min(PLHG,[],2),[],2),2);

% smoothing and normalizing
PLHGs = smoothdata(PLHG,2,'gaussian',2*Fs);
PLHGsn = 100*round((PLHGs-min(PLHGs,[],2))./max(PLHGs-min(PLHGs,[],2),[],2),2);

%% changing variable names to save time.
%PLHGsn = round(100.*(PLHGs./max(max(PLHGs))));


%% plotting the PLHG across faces on the inflated brain
% heat map showing percentage of PLHG activation.
heatMap = flipud(colormap(hot(100)));
alpha = 0.6;

% HOOD STEP
Vinds = Vinds(1:length(labels));

% % finding faces to color.
% [Fx,locX] = ismember(F(:,1),Vinds);
% [Fy,locY] = ismember(F(:,2),Vinds);
% [Fz,locZ] = ismember(F(:,3),Vinds);
%
% Finds = Fx | Fy | Fz;
%
% % what am I even doing with my life?
% A = [locX locY locZ];
% B = A(Finds,:);
% B(B==0) = NaN;
% elecs2faces = round(nanmean(B,2));

[V,F] = freesurfer_read_surf([surfPath 'lh.inflated']);


% number of nearest faces to use. 
kFaces = 250;
[Finds,D] = knnsearch(V(F),[Vinds' Vinds' Vinds'],'K',kFaces);


% converting electrodes to faces again.
PLHGfaces = zeros(0,size(PLHGsn,2));
for fc = 1:length(Vinds)
    % implements colormap smoothing
    catVar = uint8(repmat(PLHGsn(fc,:),kFaces,1).*repmat(linspace(1,0.25,kFaces)',1,size(PLHGsn,2)));
    % builds matrix of faces
    PLHGfaces = cat(1,PLHGfaces,catVar);
end
PLHGfaces(PLHGfaces==0)=1;


% movie params
framesToChangeView = 11;
origAz = 270;
origEl = 0;
propAz = 220;
propEl = -10;
moveAz = linspace(origAz,propAz,framesToChangeView);
moveEl = linspace(origEl,propEl,framesToChangeView);
nFrames = size(PLHGfaces,2);

% now making the movie.
for frm = [1:framesToChangeView framesToChangeView+1:250:nFrames]
    fprintf('\nplotting frame %d of %d...',frm,nFrames)
    % now looping over frames
    if frm==1
        % plotting surface once.whos
        [V,F,p] = plotFreeSurf_bbg(surfPath,'inflated','lh',alpha);
        
        % setting up movie stuff
        dbPath = '/media/user1/data4TB/Dropbox/Dropbox/';
        movie_file_name = [dbPath patientID 'flattenedBrain_4xspeed_V3verts.avi'];
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
        
    elseif frm<framesToChangeView %this implements a move in the view at first, so that
        % slight chagne in camera view
        view(moveAz(frm-1),moveEl(frm-1))
        
        % write frame
        writeVideo(movie_object,getframe(movie_handle));
        
        % adjust view for next section.
        view(propAz,origEl)
        
    elseif frm>=framesToChangeView
        %         try
        % selectively coloring faces.
        faceMap = zeros(size(PLHGfaces,1),3);
        for fm=1:size(PLHGfaces,1)
            faceMap(fm,:) = heatMap(PLHGfaces(fm,frm),:);
        end
        faceColor(Finds,:) = faceMap;
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


