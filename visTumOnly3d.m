function pT = visTumOnly3d(volT)
% Plots a 3D-image of a tumour (volT)

% Threshold of cell density to be considered a part of a tumour
lvlT = 0.02;

[Z,X,Y] = meshgrid(1:size(volT,2), 1:size(volT,1), 1:size(volT,3));

% Plot tumour
pT = patch(isosurface(X, Y, Z, volT, lvlT), ...
     'FaceColor', [255, 0, 0]/255, ...
     'EdgeColor','none', ...
     'FaceAlpha',1);
isonormals(volT,pT);

ax = gca;
ax.ZDir = 'reverse';
axis tight;
axis image;
axis vis3d
axis on;
box on;
grid on;
% view(-150, 15);  % [X,Y,Z] = meshgrid(...) at the beginning
view(125, 45);     % [Z,X,Y] = meshgrid(...) at the beginning 
camlight right;
camlight left;
lighting flat
rotate3d on
end
