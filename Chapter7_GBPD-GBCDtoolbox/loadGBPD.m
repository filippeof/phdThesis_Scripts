function GBPD = loadGBPD(datFname, varargin)
%Function actions:
    % Load rho, theta and mud from gbpd .dat file
%Input:
    % datFname = GBPD structure containing rho, theta and m.u.d.
%Output:
    % GBPD = Structure with coordinates and multiple of uniform distribution (m.u.d.) values 
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
%% Create gbpd structure
GBPD=dlmread(datFname);
GBPD(1,:)=[];GBPD(:,4)=[];
GBPD=array2table(GBPD);
% Rename Columns
GBPD.Properties.VariableNames{1}='theta';
GBPD.Properties.VariableNames{2}='rho';
GBPD.Properties.VariableNames{3}='mud';
% Round values
GBPD.theta=round(GBPD.theta);
% Round values and convert rho to Matlab reference frame
GBPD.rho=round(GBPD.rho-90); 

%% Area fraction calculation example
% cs=crystalSymmetry('mmm', [4.762 10.225 5.994]);
% m=Miller(vector3d('polar',GBPD.rho,GBPD.theta),cs);% transform polar coordinates to miller idx
% mr=round(m); %round miller index 
% 
% mt=table(m.hkl);% create table from hkl  indexes...
% mt.Properties.VariableNames{1}='hkl'; %rename column
%  mrt=table(mr.hkl);
%  mrt.Properties.VariableNames{1}='hkl_Round';
% 
% areafraction=(GBPD.mud/81)*100; 
% dists=(areafraction.^-1)/max((areafraction.^-1)); %distance normalized from 0 to 1
% 
% GBPD=[GBPD table(areafraction) table(dists) mt mrt];
%%
fprintf('\n -loadGBPD: GBPD imported with the following properties: \n')
disp(GBPD.Properties.VariableNames)
end

