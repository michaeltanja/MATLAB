% Michael Tanja and Jonathan DiBacco

% The first segment of the code is where all the parameters are stated. 
% The only thing you will have to change is the NperRow and the NperCol
% The 'dx' changes according to the inputs for NperRow/Col


clear all
close all

%%
%parameters and
%mesh properties
prompt1 = 'What resolution of dx would do you want? ';
prompt2 = 'What is the length and width of this square plate? ';
L  = input(prompt2);         %unsure about what this does, but it was in the example
W  = L;         %unsure about what this does, but it was in the example
dx = input(prompt1);
qx = 100000;    %heat flux bottom
qy = -65000;     %heat flux side
NperRow = (L/dx) + 1; % Nodes per row
NperCol = NperRow; % Nodes per column
nNodes = NperRow * NperCol;
A = zeros(NperRow, NperCol);
b = zeros(NperRow,1);

k  = 250;       %conduction coefficient
h  = 250;       %convection coefficient
Tinf = 0;       %fluid temp in celsius


%next in line is node set up
%% 
% Center Node/s

for ii = NperRow+2 : ((2*NperRow-1))
    for jj = ii:NperRow:nNodes-(NperRow+1)%NperRow+2 : (NperRow^2) - (NperRow+1) %put this first in order, matlab reads top down
    
    A(jj, jj+NperRow) = 1; %upper
    A(jj, jj-NperRow) = 1; %lower
    A(jj, jj+1) = 1;       %right
    A(jj, jj-1) = 1;       %left
    A(jj, jj) = -4;        %center
    b(jj) = 0;             %boundary
    end
end

%%
% Convection node/s on the right

for ii = 2*NperRow : NperRow : (NperRow^2) - NperRow              %needs to be more generalized
    jj=ii
    A(ii, jj+NperRow) = 1;       %upper
    A(ii, jj-NperRow) = 1;       %lower
                                 %no right
    A(ii, jj-1) = 2;             %left
    A(ii, jj) = (-4-((2*h*dx)/k)); %center
    b(ii) = -((2*h*dx)/k)*Tinf;   %boundary
    
end

%% 
% Convection node/s on top
for ii = (nNodes )-(NperRow - 2) : (nNodes -1) %generalized
    jj= ii
                                                          %no top term
    A(ii, jj-NperRow) = 2;                                %lower
    A(ii, jj+1) = 1;                                      %right 
    A(ii, jj-1) = 1;                                      %left
    A(ii, jj) = (-4-((2*h*dx)/k));                    %center
    b(ii) = -(2*h*dx*Tinf)/k;                              %boundary
end

%starting corners 1-4
%%
% Corner 1 (this boundary condition should have same initiation
%   for any matrix size)
for ii = (NperCol*NperRow)-(NperRow-1)     %generalized
    jj=ii
                                           %no top term 
    A(ii, jj-NperRow) = 1;                 %lower
    A(ii, jj+1) = 1;                       %right 
                                           %no left term
    A(ii,jj) = (-2-((2*h*dx)/k));       %center
    b(ii) = -(((2*h*dx*Tinf)/k) +((2*qy*dx)/k)); %boundary
                       %EDIT: change flow of heat on boundary condition
end

%%
% Corner 2 (This shouldn't change either)

for ii = nNodes
    jj=ii
                                     %no top term
    A(ii, jj-NperRow) = 1;           %lower term
                                     %no right term
    A(ii, jj-1) = 1;                 %left
    A(ii, jj) = (-(2*h*dx)/(k)) - 2; %center
    b(ii) = -(2*h*dx*Tinf)/k;              %boundary
end

%%
% Corner 3 (should stay the same for all meshes)

for ii = 1
    jj=ii
    A(ii, jj+NperRow) = 1;  %upper
                            %no lower
    A(ii, jj+1) = 1;        %right
                            %no left term
    A(ii, jj) = -2;         %center
    b(ii) = -((dx/k)*(qy+qx)); %boundary
               %EDIT: Changed flow of heat on the boundary condition
end

%%
% Corner 4 (should be the same for all matrices)

for ii = NperRow
    jj=ii
    A(ii, jj+NperRow) = 1;                 %upper
                                           %no lower term
                                           %no right term
    A(ii, jj-1) = 1;                       %left term
    A(ii,jj) = -2-((h*dx)/k);              %center
    b(ii) = -(((qx*dx) + (h*dx*Tinf))/k);  %boundary
end

%corners are finished
%Start heat sources
%%
%heat source left side

for ii = NperRow+1 : NperRow : (nNodes)-(NperRow+1)                 %try to generalize for all matrices
    jj=ii
    A(ii, jj+NperRow) = 1; %upper
    A(ii, jj-NperRow) = 1; %lower
    A(ii, jj+1) = 2;       %right
                           %no left term
    A(ii, jj) = -4;        %center
    b(ii) = -(2*(qy)*dx)/k;   %boundary
             %EDIT: Changed flow of heat on the boundary condition, for qy
end

%%
%heat source bottom

for ii = 2: (NperRow-1)    %generalized
    jj=ii
    A(ii, jj+NperRow) = 2; %upper
                           %no lower term
    A(ii, jj+1) = 1;       %right term
    A(ii, jj-1) = 1;       %left term
    A(ii, jj) = -4;        %center
    b(ii) = -(2*qx*dx)/k;   %boundary
end

%%
%Matrices are all set up, should be able to solve and graph now
%v = b'; %turned array into vector to make dimensions agree
T = A\b;

%%
%creating the mesh plot 

x = [0: dx :L]; %x-axis grid
y = [0: dx :W]; %y-axis grid

[X,Y] = meshgrid(x,y);
%Create a loop that converts the vector back into a mesh to plot
%takes the temperature vector and reorganizes into a matrix

nIndex = 1; %this index will count down the temperature
for aa = 1:NperRow %nNodes/NperRow
    for bb = 1:NperRow
        TMesh(aa,bb) = T(nIndex);
        nIndex = nIndex +1;
    end
end
%TMesh1 = abs(TMesh);
%recycled some of the example code
z = length(x);
contourf(X,Y,TMesh)
%plot(x,TMesh((z+1)/2, :))
%plot(y,TMesh(:,(z+1)/2))
colorbar
title('Contour plot of temperature distribution(C) as a function of position(m)')
%contour looks ok,but heat might be flowing in the wring direction
%change the sign for heat flow on the right face of the figure
%^this should affect values for T1, T4, T7
