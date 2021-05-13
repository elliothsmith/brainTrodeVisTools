% Brush
load CentroidCoordinates.mat

% Create X,Y,Z, variable for ease of reading code below
X = CentroidCoordinates(:,1);
Y = CentroidCoordinates(:,2);
Z = CentroidCoordinates(:,3);

% Create three 2-D subplo ts
H = figure

    subplot(1,3,1), scatter(X,Y)
    subplot(1,3,2), scatter(X,Z)
    subplot(1,3,3), scatter(Y,Z)
    linkdata
    
    brush on
    brushobj = brush(H)

    