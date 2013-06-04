function nR = meanRes(R,ncellY,ncellX,nlayer)
nR = zeros(ncellX,ncellY,nlayer);
for l = 1:nlayer
   % fprintf('calculating cell centered resistivity of layer %i...     \n',l)
    for j = 1:ncellY
        for i = 1:ncellX
            nR(i,j,l) = mean ([R(i,j,l),R(i+1,j,l),R(i+1,j+1,l),R(i,j+1,l)]);
        end
    end
end
end
