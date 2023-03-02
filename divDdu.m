function diffOp = divDdu(u, Davg, OmegaEdges,p)

if ~exist('p', 'var')
    p=2;
end

udifPos = [1;-1;0];
udifNeg = [0;1;-1];

% Direction 1 (y)
% u1Pos = u(x, y+1, z) - u(x, y, z), and = 0 if (x,y+1,z) is outside Omega
% u1Neg = u(x, y, z) - u(x, y-1, z), and = 0 if (x,y-1,z) is outside Omega
u1Pos = convn(u, reshape(udifPos, 3, 1, 1), 'same') .* OmegaEdges.Pos1;
u1Neg = convn(u, reshape(udifNeg, 3, 1, 1), 'same') .* OmegaEdges.Neg1;

% Direction 2 (x)
% u2Pos = u(x+1, y, z) - u(x, y, z), and = 0 if (x+1,y,z) is outside Omega
% u2Neg = u(x, y, z) - u(x-1, y, z), and = 0 if (x-1,y,z) is outside Omega
u2Pos = convn(u, reshape(udifPos, 1, 3, 1), 'same') .* OmegaEdges.Pos2;
u2Neg = convn(u, reshape(udifNeg, 1, 3, 1), 'same') .* OmegaEdges.Neg2;

% Direction 3 (z)
% u3Pos = u(x, y, z+1) - u(x, y, z), and = 0 if (x,y,z+1) is outside Omega
% u3Neg = u(x, y, z) - u(x, y, z-1), and = 0 if (x,y,z-1) is outside Omega
u3Pos = convn(u, reshape(udifPos, 1, 1, 3), 'same') .* OmegaEdges.Pos3;
u3Neg = convn(u, reshape(udifNeg, 1, 1, 3), 'same') .* OmegaEdges.Neg3;

duPosNorm = u1Pos.^2 + u2Pos.^2 + u3Pos.^2;
duNegNorm = u1Neg.^2 + u2Neg.^2 + u3Neg.^2;

duPosNorm(duPosNorm==0)=1;
duNegNorm(duNegNorm==0)=1;


diffOp = Davg.Pos1 .* u1Pos.*duPosNorm.^(p/2-1) - Davg.Neg1 .* u1Neg.*duNegNorm.^(p/2-1)  + ...
         Davg.Pos2 .* u2Pos.*duPosNorm.^(p/2-1)  - Davg.Neg2 .* u2Neg.*duNegNorm.^(p/2-1)  + ...
         Davg.Pos3 .* u3Pos.*duPosNorm.^(p/2-1)  - Davg.Neg3 .* u3Neg.*duNegNorm.^(p/2-1) ;
end
