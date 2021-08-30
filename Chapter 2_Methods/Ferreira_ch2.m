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
%% Chapter 2: methods

% Pressure and Temperature range of experiments
close all
colors=lines;
r = randperm(size(colors,1)); % permute row numbers
colors = colors(r,:);
labels_app={'Patterson Apparattus', 'Piston Cylinder', 'MAVO press','Multi-Anvil'};
P_range= [0 0.5; 0 4; 0 25; 0 40];
T_range= [0 1300; 0 1800; 0 2000; 0 2100];

for i=length(labels_app):-1:1
    rectangle('Position',[P_range(i,1) T_range(i,1) P_range(i,2) T_range(i,2)],'FaceColor',[colors(i,:),0.2])
    if i==1 || i==2
        text(mean(P_range(i,1:2)), mean(T_range(i,1:2)),labels_app{i},'Rotation',90)
    elseif i==3
        text(mean(P_range(i,1:2)), mean(T_range(i,1:2)),labels_app{i})
    else
        text(mean(P_range(i,1:2))+mean(P_range(i-1,1:2)), mean(T_range(i,1:2)),labels_app{i})
    end
end

axis tight
xlabel('Pressure (GPa)','FontSize',14); ylabel('Temperature (°C)','FontSize',14);
%
h=[40:10: 410]; %depth
P=0.0357.*h-0.61;%Pressure(GPa)
T=(-4.9960e-07)*(h.^4)+(5.7589e-04)*(h.^3)-0.24*(h.^2)+42.51.*h-1181.58;
hold on
p=plot(P,T,'-.','LineWidth',4,'Color','k');
p.Color(4) = 0.5;
hold off
print('/Users/fereira/Documents/DOC/GBCD_Olivine/Thesis/Illustrations/PxT_Earth','-dpdf', '-painters')
%%
% Olivine seismic anisotropy
cs_tensor = crystalSymmetry('mmm',[4.7646,10.2296,5.9942],...
  'x||a','z||c','mineral','Olivine');

M = [[320.5  68.15  71.6     0     0     0];...
    [ 68.15  196.5  76.8     0     0     0];...
    [  71.6   76.8 233.5     0     0     0];...
    [   0      0      0     64     0     0];...
    [   0      0      0      0    77     0];...
    [   0      0      0      0     0  78.7]];
rho=3.355;
C = stiffnessTensor(M,cs_tensor,'density',rho);
[vp,vs1,vs2,pp,ps1,ps2] = C.velocity('harmonic');

% plotting convention - plot a-axis to east
plota2east;

% set colour map to seismic color map : blue2redColorMap
setMTEXpref('defaultColorMap',blue2redColorMap)

% some options
blackMarker = {'Marker','s','MarkerSize',10,'antipodal',...
  'MarkerEdgeColor','white','MarkerFaceColor','black','doNotDraw'};
whiteMarker = {'Marker','o','MarkerSize',10,'antipodal',...
  'MarkerEdgeColor','black','MarkerFaceColor','white','doNotDraw'};

% some global options for the titles
%titleOpt = {'FontSize',getMTEXpref('FontSize'),'visible','on'}; %{'FontSize',15};
titleOpt = {'visible','on','color','k'};

% Setup multiplot
% define plot size [origin X,Y,Width,Height]
mtexFig = mtexFigure('position',[0 0 1000 1000]);
%**************************************************************************
% Vp : Plot P-wave velocity (km/s)
%**************************************************************************

% Plot P-wave velocity (km/s)
plot(vp,'contourf','complete','upper')
mtexTitle('Vp (km/s)',titleOpt{:})

% extrema
[maxVp, maxVpPos] = max(vp);
[minVp, minVpPos] = min(vp);

% percentage anisotropy
AVp = 200*(maxVp-minVp) / (maxVp+minVp);

% mark maximum with black square and minimum with white circle
hold on
plot(maxVpPos.symmetrise,blackMarker{:})
plot(minVpPos.symmetrise,whiteMarker{:})
hold off

% subTitle
xlabel(['Vp Anisotropy = ',num2str(AVp,'%6.1f')],titleOpt{:})

%%
cs = crystalSymmetry('mmm',[4.7646,10.2296,5.9942],...
  'x||a','z||c','mineral','Olivine');

h = Miller({1,0,0},{0,1,0},{0,0,1},cs,'hkl'); %Poles
v=vector3d(h)./norm(vector3d(h));
vp=velocity(C,v,rho);

csP=crystalShape.olivine
csP.N
csP.plot