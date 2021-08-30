function  GBCD_par=getGBCDpar(ebsd, cs)
%Function actions:
    %Reads EBSD file and outputs:
    % msym, rot, rad, resolution, seg, stepSize, cutOff
    % to input.txt files. Edit this file for the desired parameters.
%Input:
    % ebsd = root folder with gbpd/gbcd programs
    % segments_fname = file name of segments
    % GBCD_par = gbpd/gbcd/parameters (obtained from getGBCDpar.m)
%Output:
    %datFname: Path to the GBPD .dat file located in the graph_pd folder to be ploted
%{
    Writen by Filippe Ferreira (BGI, Uni-Bayreuth, Germany) supervised by Katharina Marquardt (Imperial College London) and Marcel
    Thielmann (BGI, Uni-Bayreuth, Germany). 

    Copyright 2021 Ferreira & Marquardt
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the
    Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
    and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%}
%% Check symmetry
switch cs.lattice
    case 'cubic'
        msym=3;
    case 'hexagonal'
        msym=2;
    case 'tetragonal'
        msym=1;
    case 'trigonal'
        msym=4;
    case 'orthorhombic'
        msym=5;
    otherwise
        error('Crystal symmetry not available')
end

% Rotation %** Change accordingly
%  rot=0;% HKL system
rot=1; % TSL system
%TODO: identify system?

% Euler angles in rad or deg 
% rad=0% degree
rad=1; % rad

% Resolution (bins per 90 degrees) -maximum:18.     %** Change accordingly
resolution=[9, 9];

% Segment analysys %** Change accordingly
% Recomended to leave at zero when using Mtex exported data
% (in mtex every segment has the same size)
seg= 0;

%Check step size
ucell=ebsd.unitCell;
area_shape=polyarea(ucell(:,1),ucell(:,2));

switch size(ucell,1)
    case 4 % square grid
        stepSize=sqrt(area_shape); 
    case 6 % hexagonal grid
        stepSize=sqrt(area_shape/sind(60));
    case 16 % ? grid
        stepSize=sqrt((area_shape)/pi);% check this
    otherwise
         stepSize= input('Enter EBSD step size and press enter:  ');
end

% small length cutoff %** Change accordingly
% Recomended to leave at zero when using Mtex exported data
cutOff= 0.0;

%Join all parameters in array 
GBCD_par=[msym, rot, rad, resolution, seg, stepSize, cutOff];
%
fprintf('\n -getGBCDpar: The following parameters will be used in the GBPD/CD calculations: \nSymmetry | Rotation | Radians | Resolution | Segment analysys | Step size | Small length Cut off \n  %11.0f | %12.0f | %11.0f | %10.0f x %.0f| %25.0f | %10.2f | %28.2f  \n ', GBCD_par);
