%Create colormap.dat from Paraview preset colorbar
%To include Colormap
%1. Open Paraview Edit Color Map
%2. Export with .json
%3. Run this function to create colormap.dat
%4. OKUBO 05/2016
function TranslateCMap
    delete colormap.dat
    A = importdata('BuRd.json',',');

    for i = 1:(length(A.data)/4)
        temp(i,1) = A.data(4*(i-1)+1);
        temp(i,2) = A.data(4*(i-1)+2)*255;
        temp(i,3) = A.data(4*(i-1)+3)*255;
        temp(i,4) = A.data(4*(i-1)+4)*255;
    end

    columnum = length(A.data)/4;

    save('colormap.dat','columnum','temp','-ascii');
end