function pTnew = visTum3dUpdate(volT, pTold)
% Updates a 3D-image of brain matter, by removing previous tumour and
% plotting the a new Tumor volume instead

lvlT = 0.16;

if exist('pTold', 'var')
    delete(pTold);
end

% Plot tumour
hold on;

[Z,X,Y] = meshgrid(1:size(volT,2), 1:size(volT,1), 1:size(volT,3));

pTnew = patch(isosurface(X, Y, Z, volT, lvlT), ...
             'FaceColor', [255, 0, 0]/255, ...
             'EdgeColor','none', ...
             'FaceAlpha',1);
isonormals(volT, pTnew);
hold off
end
