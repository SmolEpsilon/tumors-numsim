function pT = visTum3d(volD, volT, cutCoord)
% Plots a 3D-image of brain matter (volD) and a tumour (volT)
% where the brain matter cut off along 3rd dimension at coordinates
% given by cutCoord

clf; 

lvlGW = 0.1;
lvlT = 0.02;

[Z,X,Y] = meshgrid(1:size(volD,2), 1:size(volD,1), 1:size(volD,3));
Xt = X; Xt(volD<lvlGW) = NaN;
Yt = Y; Yt(volD<lvlGW) = NaN;
Zt = Z; Zt(volD<lvlGW) = NaN;

% prepare caps
Xcap1 = squeeze(Xt(cutCoord(1),1:cutCoord(2),:));
Ycap1 = squeeze(Yt(cutCoord(1),1:cutCoord(2),:));
Zcap1 = squeeze(Zt(cutCoord(1),1:cutCoord(2),:));
Gcap1 = squeeze(volD(cutCoord(1),1:cutCoord(2),:));

Xcap2 = squeeze(Xt(cutCoord(1):end,cutCoord(2),:));
Ycap2 = squeeze(Yt(cutCoord(1):end,cutCoord(2),:));
Zcap2 = squeeze(Zt(cutCoord(1):end,cutCoord(2),:));
Gcap2 = squeeze(volD(cutCoord(1):end,cutCoord(2),:));

% prepare coordinates for cut
X(cutCoord(1)+1:end, 1:cutCoord(2)-1, :) = NaN;
Y(cutCoord(1)+1:end, 1:cutCoord(2)-1, :) = NaN;
Z(cutCoord(1)+1:end, 1:cutCoord(2)-1, :) = NaN;

[Zt,Xt,Yt] = meshgrid(1:size(volD,2), 1:size(volD,1), 1:size(volD,3));


% Plot gray+white matter
% pGW = ...
patch(isosurface(X,Y,Z,volD, lvlGW), ...
     'FaceColor', [255,173,96]/255, ...
     'EdgeColor','none', ...
     'FaceAlpha',0.6);
% isonormals(volD,pGW);
hold on

% Cap the cut region
surf(Xcap1, Ycap1, Zcap1, Gcap1, 'FaceColor', 'interp', 'EdgeColor', 'none');
surf(Xcap2, Ycap2, Zcap2, Gcap2, 'FaceColor', 'interp', 'EdgeColor', 'none');

% Plot tumour
pT = patch(isosurface(Xt,Yt,Zt,volT, lvlT), ...
     'FaceColor', [255, 0, 0]/255, ...
     'EdgeColor','none', ...
     'FaceAlpha',1);
isonormals(volT,pT);

hold off
colormap gray;
ax = gca;
ax.CLim = [-max(volD(:)), max(volD(:))];
ax.ZDir = 'reverse';
axis tight;
axis image;
axis vis3d
axis off;
% view(-150, 15);  % if [X,Y,Z] = meshgrid(...) is at the beginning
view(125, 45);     % if [Z,X,Y] = meshgrid(...) is at the beginning 
camlight right;
camlight left;
lighting gouraud
rotate3d on
