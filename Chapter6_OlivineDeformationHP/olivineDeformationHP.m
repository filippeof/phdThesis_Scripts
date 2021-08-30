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
%% Olivine deformation at HP


% Sample HH221 (FSGHPD1A)
    fname_HH221_SC= '/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH221/HH221-SC/map20190511174102326_crop_clean.ang';
    fname_HH221_SG= '/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH221/HH221-SG/map20190510212038753_crop_clean.ang';
 
    % Cleaning parameters
    % Minimum confidence index
    minCI =0.1; 
    % Minimum  grain size (points)
    minGS =20; 
    % Minimum grain misorientation angle threshold
    minANG =20*degree; 
    % Grain boundaries smooth factor (number of iterations)
    sF=5;
    %ODF halfwidth
    hW=8*degree;
    % EBSD phase name chosen for plots
    phaseName='olivine'; 
    %Set phase colors
    ol_color=[97, 148, 63]./255;
    en_color=[98, 63, 35]./255;

    % Set coordinate reference frame
    setMTEXpref('xAxisDirection','east');setMTEXpref('zAxisDirection','intoPlane'); %EDAX: A1 up(-y), A2 left(-x)
    cs=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');
    oM = ipfHSVKey(cs);
    oM.inversePoleFigureDirection = zvector;% 001 IPF colorcoding 

    ebsdFileList= {fname_HH221_SC, fname_HH221_SG};
    % Loop over all ebsd data, save grains data and plot maps      
    for i=1:length(ebsdFileList)
        ebsdFname=ebsdFileList{i};
        [fPt,sampleName,~] = fileparts(ebsdFname);
        grains_fname=fullfile(fPt,[sampleName,'_grains.mat']);
        ebsd_fname=fullfile(fPt,[sampleName,'_ebsdClean.mat']);

        % ------Import data ---------
        % Crystal symmetries of all phases
        CS=getCS(ebsdFname);
        % Import EBSD
        ebsd_raw =EBSD.load(ebsdFname,'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');
        ebsd=ebsd_raw;

        % ------Clean data-----------
        % Remove data with CI lower than the defined minimum
        try
            idNotIdx=unique(ebsd('notIndexed').phase);
        catch 
        	idNotIdx=0;
        end

        ebsd(ebsd.prop.ci<minCI).phase=idNotIdx;
        % Calculate grains
        [grains,ebsd.grainId] = calcGrains(ebsd,'threshold',minANG);
        %Remove holes
        try
            notIndexed = grains('notIndexed');
            toRemove = notIndexed(notIndexed.grainSize ./ notIndexed.boundarySize<0.8);
            ebsd(toRemove) = [];
        catch
        end
        % Remove ebsd data of grains smaller than the defined minimum grain size
        ebsd(grains(grains.grainSize<minGS)) = [];
        % Recalculate grains after cleaning 
        [grains,ebsd.grainId] = calcGrains(ebsd,'threshold',minANG);
%          Smooth grain boundaries
        grains=smooth(grains,sF);

%         % Remove outer boundary grains
%         outerBoundary_id = any(grains.boundary.grainId==0,2);
%         grain_id = grains.boundary(outerBoundary_id).grainId;
%         grain_id(grain_id==0) = [];
%         grains('id', grain_id) = [];

        % ------ Plot data----------
        %IPF
        figure(99);
        color = oM.orientation2color(ebsd(phaseName).orientations);
        setMTEXpref('figSize','large');
        setMTEXpref('outerPlotSpacing',0);
        plot(ebsd_raw,ebsd_raw.prop.iq,'micronbar','on')
        mtexColorMap black2white
        hold on; plot(ebsd(phaseName), color, 'FaceAlpha',0.85);% plot orientations
        hold on; plot(grains('indexed').boundary,'k', 'linewidth', 0.9, 'edgealpha', 0.6);hold off;
        saveFigure(fullfile(fPt, [sampleName,'_ipf.png']))

        %Grains
        figure(100);
        setMTEXpref('figSize','large');
        setMTEXpref('outerPlotSpacing',0);
        plot(ebsd_raw,ebsd_raw.prop.iq,'micronbar','on')
        mtexColorMap gray
        grains('ol').color=ol_color; 
        if ~isempty(grains('en')); grains('en').color=en_color; end
        
        hold on; plot(grains('indexed'), 'FaceAlpha',0.85);% plot orientations
        hold on; plot(grains('indexed').boundary,'k', 'linewidth', 0.9, 'edgealpha', 0.6);hold off;
      saveFigure(fullfile(fPt, [sampleName,'_phaseMap.png']))

        % ------- Save data-----------
%         save (ebsd_fname, 'ebsd');
%         save (grains_fname, 'grains');
    end
    %%  Load data HH221
    ebsd_SC=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH221/HH221-SC/map20190511174102326_crop_clean_ebsdClean.mat');
    ebsd_SC=ebsd_SC.ebsd;

    ebsd_SG=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH221/HH221-SG/map20190510212038753_crop_clean_ebsdClean.mat');
    ebsd_SG=ebsd_SG.ebsd;
    
    grains_SC=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH221/HH221-SC/map20190511174102326_crop_clean_grains.mat');
    grains_SC=grains_SC.grains;

    grains_SG=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH221/HH221-SG/map20190510212038753_crop_clean_grains.mat');
    grains_SG=grains_SG.grains;


%%
figure;
ipfKey = axisAngleColorKey(ebsd_SC('ol'));
% set the grain mean orientations as reference orinetations
ipfKey.oriRef = grains_SC(ebsd_SC('ol').grainId).meanOrientation;
% colorKey.oriRef = grains_SC(ebsd_SC.grainId).meanOrientation;

% lets plot the result
plot(ebsd_SC('ol'),ipfKey.orientation2color(ebsd_SC('ol').orientations))
hold on
plot(grains_SC.boundary,'linewidth',2)
hold off
% orientation plot
% ipfKey.oriRef = grains_SC(ebsd_SC_S('ol').grainId).meanOrientation;
% plot(ebsd_SC_S('ol'),ebsd_SC_S('ol').orientations)
% hold on
% plot(grains_SC.boundary,'linewidth',2)
% hold off
% saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH221_IPF_SC','-m2.5')
%% KAM-SC

close all
F = halfQuadraticFilter;
% smooth the data
ebsd_SC_S = smooth(ebsd_SC,F,'fill',grains_SC);
setMTEXpref('outerPlotSpacing',10);
setMTEXpref('innerPlotSpacing',10);
kam = KAM(ebsd_SC_S('ol'),'order',2,'threshold',5*degree);
figure
f=newMtexFigure('figsize','tiny');
plot(ebsd_SC_S('ol'), kam./degree,'micronbar','on');
mtexColorbar('title', 'KAM [º]')
hold on
plot(grains_SC('ol').boundary,'linewidth',1.5,'micronbar','on')
hold off
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH221_KAM_SC')

%% KAM-SG
f=newMtexFigure('figsize','tiny');;
setMTEXpref('outerPlotSpacing',10);
setMTEXpref('innerPlotSpacing',10);
F = halfQuadraticFilter;
% smooth the data
ebsd_SG_S = smooth(ebsd_SG,F, 'fill', grains_SG); 

kam = KAM(ebsd_SG_S('ol'),'order',2,'threshold',5*degree);
figure
plot(ebsd_SG_S('ol'), kam./degree,'micronbar','on');
mtexColorbar('title', 'KAM [º]')
hold on
plot(grains_SG('ol').boundary,'linewidth',1.5,'micronbar','on')
hold off
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH221_KAM_SG')

%% mis2mean-SC
[grainsSC_S,ebsd_SC_S.grainId,ebsd_SC_S.mis2mean] = calcGrains(ebsd_SC_S,'threshold',20*degree);
%%
close all;figure;
f=newMtexFigure('figsize','tiny');
[~,mP] =plot(ebsd_SC_S('ol'), ebsd_SC_S('ol').mis2mean.angle./degree,'micronbar','on');
CLim(gcm, [0 10])
mtexColorbar('title', 'Misorientatation to mean orientation [º]')
hold on
plot(grainsSC_S('ol').boundary,'linewidth',1.5,'micronbar','on')
hold off
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH221_mis2mean_SC')
mtexColorbar('off')
mP.micronBar.visible = 'off';
% 3 insets
xtd=ebsd_SC_S.extend;

xmax=xtd(2);
xmid=round((xtd(2)-xtd(1))/7);
ymid=round((xtd(4)-xtd(3))/2);
ivl=125;
axis tight
for i=[2,4,6]
    pos=[i*xmid-ivl*2 ymid-ivl ivl*2 ivl*2];
    rectangle('Position',pos,'EdgeColor','w')
    axis([pos(1) pos(1)+ivl*2 pos(2) pos(2)+ivl*2])
%     az=180; el=-90; view(az, el);
    set(gca,'Ydir','reverse')
    CLim(gcm, [0 10])
    print(gcf, '-dpng', '-r300', ['/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH221_m2mean_SC_inset_',num2str(i),'.png']);
end

% saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH221_mis2mean_SC')

%% mis2mean-SG
[grainsSG_S,ebsd_SG_S.grainId,ebsd_SG_S.mis2mean] = calcGrains(ebsd_SG_S,'threshold',20*degree);
close all;figure;
f=newMtexFigure('figsize','tiny'); 
plot(ebsd_SG_S('ol'), ebsd_SG_S('ol').mis2mean.angle./degree,'micronbar','on')
mtexColorbar('title', 'Misorientatation to mean orientation [º]')
hold on
plot(grainsSG_S('ol').boundary,'linewidth',1.5,'micronbar','on')
hold off
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH221_mis2mean_SG')
%%  ____________________
%HH222
%
    fname_HH222_SC= '/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH222/HH222-SC/map20200825174846366_Rescan_clean_ps.ang';
    fname_HH222_SG= '/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH222/HH222-SG/map20200826214013815_Rescan_clean_ps.ang';
 
  % Cleaning parameters
    % Minimum confidence index
    minCI =0.1; 
    % Minimum  grain size (points)
    minGS =20; 
    % Minimum grain misorientation angle threshold
    minANG =20*degree; 
    % Grain boundaries smooth factor (number of iterations)
    sF=5;
    %ODF halfwidth
    hW=8*degree;
    % EBSD phase name chosen for plots
    phaseName='olivine'; 
    %Set phase colors
    ol_color=[97, 148, 63]./255;
    en_color=[98, 63, 35]./255;

    % Set coordinate reference frame
    setMTEXpref('xAxisDirection','east');setMTEXpref('zAxisDirection','intoPlane'); %EDAX: A1 up(-y), A2 left(-x)
    cs=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');
    oM = ipfHSVKey(cs);
    oM.inversePoleFigureDirection = zvector;% 001 IPF colorcoding 

    ebsdFileList= {fname_HH222_SC, fname_HH222_SG};
    % Loop over all ebsd data, save grains data and plot maps      
    for i=1:length(ebsdFileList)
        ebsdFname=ebsdFileList{i};
        [fPt,sampleName,~] = fileparts(ebsdFname);
        grains_fname=fullfile(fPt,[sampleName,'_grains.mat']);
        ebsd_fname=fullfile(fPt,[sampleName,'_ebsdClean.mat']);

        % ------Import data ---------
        % Crystal symmetries of all phases
        CS=getCS(ebsdFname);
        % Import EBSD
        ebsd_raw =EBSD.load(ebsdFname,'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');
        ebsd=ebsd_raw;

        % ------Clean data-----------
        % Remove data with CI lower than the defined minimum
        try
            idNotIdx=unique(ebsd('notIndexed').phase);
        catch 
        	idNotIdx=0;
        end

        ebsd(ebsd.prop.ci<minCI).phase=idNotIdx;
        % Calculate grains
        [grains,ebsd.grainId] = calcGrains(ebsd,'threshold',minANG);
        %Remove holes
        try
            notIndexed = grains('notIndexed');
            toRemove = notIndexed(notIndexed.grainSize ./ notIndexed.boundarySize<0.8);
            ebsd(toRemove) = [];
        catch
        end
        % Remove ebsd data of grains smaller than the defined minimum grain size
        ebsd(grains(grains.grainSize<minGS)) = [];
        % Recalculate grains after cleaning 
        [grains,ebsd.grainId] = calcGrains(ebsd,'threshold',minANG);
%          Smooth grain boundaries
        grains=smooth(grains,sF);
        
       %         % Remove outer boundary grains
%         outerBoundary_id = any(grains.boundary.grainId==0,2);
%         grain_id = grains.boundary(outerBoundary_id).grainId;
%         grain_id(grain_id==0) = [];
%         grains('id', grain_id) = [];

        % ------ Plot data----------
        %IPF
        figure(99);
        color = oM.orientation2color(ebsd(phaseName).orientations);
        setMTEXpref('figSize','large');
        setMTEXpref('outerPlotSpacing',0);
        plot(ebsd_raw,ebsd_raw.prop.iq,'micronbar','on')
        mtexColorMap black2white
        hold on; plot(ebsd(phaseName), color, 'FaceAlpha',0.85);% plot orientations
        hold on; plot(grains('indexed').boundary,'k', 'linewidth', 0.9, 'edgealpha', 0.6);hold off;
        saveFigure(fullfile(fPt, [sampleName,'_ipf.png']))

        %Grains
        figure(100);
        setMTEXpref('figSize','large');
        setMTEXpref('outerPlotSpacing',0);
        plot(ebsd_raw,ebsd_raw.prop.iq,'micronbar','on')
        mtexColorMap gray
        grains('ol').color=ol_color; 
        if ~isempty(grains('en')); grains('en').color=en_color; end
        
        hold on; plot(grains('indexed'), 'FaceAlpha',0.85);% plot orientations
        hold on; plot(grains('indexed').boundary,'k', 'linewidth', 0.9, 'edgealpha', 0.6);hold off;
        saveFigure(fullfile(fPt, [sampleName,'_phaseMap.png']))

        % ------- Save data-----------
%         save (ebsd_fname, 'ebsd');
%         save (grains_fname, 'grains');
    end
    %%
    %HH222 Load data
    ebsd_SC=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH222/HH222-SC/map20200825174846366_Rescan_clean_ps_ebsdClean.mat');
    ebsd_SC=ebsd_SC.ebsd;

    ebsd_SG=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH222/HH222-SG/map20200826214013815_Rescan_clean_ps_ebsdClean.mat');
    ebsd_SG=ebsd_SG.ebsd;
    
    grains_SC=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH222/HH222-SC/map20200825174846366_Rescan_clean_ps_grains.mat');
    grains_SC=grains_SC.grains;

    grains_SG=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH222/HH222-SG/map20200826214013815_Rescan_clean_ps_grains.mat');
    grains_SG=grains_SG.grains;

%% KAM-SC
F = halfQuadraticFilter;
% smooth the data
ebsd_SC_S = smooth(ebsd_SC,F,'fill',grains_SC); 

kam = KAM(ebsd_SC_S('ol'),'order',2,'threshold',4*degree);
figure
f=newMtexFigure('figsize','tiny'); 

plot(ebsd_SC_S('ol'), kam./degree);
mtexColorbar
hold on
plot(grains_SC('ol').boundary,'linewidth',2)
hold off
%% KAM-SG
figure
f=newMtexFigure('figsize','tiny'); 
setMTEXpref('outerPlotSpacing',10);
setMTEXpref('innerPlotSpacing',10);
F = halfQuadraticFilter;
% smooth the data
ebsd_SG_S = smooth(ebsd_SG,F, 'fill', grains_SG); 
ipfKey.oriRef = grains_SG(ebsd_SG_S('ol').grainId).meanOrientation;

kam = KAM(ebsd_SG_S('ol'),'order',2,'threshold',4*degree);
figure
plot(ebsd_SG_S('ol'), kam./degree);
mtexColorbar
hold on
plot(grains_SG('ol').boundary,'linewidth',2)
hold off
%% mis2mean-SC
[grainsSC_S,ebsd_SC_S.grainId,ebsd_SC_S.mis2mean] = calcGrains(ebsd_SC_S,'threshold',20*degree);
%%
close all;figure;
f=newMtexFigure('figsize','tiny');
[~,mP] =plot(ebsd_SC_S('ol'), ebsd_SC_S('ol').mis2mean.angle./degree,'micronbar','on');
mP.micronBar.length=150;
CLim(gcm, [0 10])
mtexColorbar('title', 'Misorientatation to mean orientation [º]')
hold on
plot(grainsSC_S('ol').boundary,'linewidth',1.5,'micronbar','on')
hold off
% 3 insets
xtd=ebsd_SC_S.extend;
xmid=round((xtd(2)-xtd(1))/7);
ymid=round((xtd(4)-xtd(3))/2);
sqSZ=150;

for i=[2,4,6]
    pos=[i*xmid-sqSZ ymid-sqSZ/2 sqSZ sqSZ];
    rectangle('Position',pos,'EdgeColor',[241 89 43]/255,'LineStyle','--','LineWidth',2)
end
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH222_mis2mean_SC')
mtexColorbar('off')
mP.micronBar.visible = 'off';
axis square tight
for i=[2,4,6]
    pos=[i*xmid-sqSZ ymid-sqSZ/2 sqSZ sqSZ];
    axis([pos(1) pos(1)+sqSZ pos(2) pos(2)+sqSZ])
%     az=180; el=-90; view(az, el);
    set(gca,'Ydir','reverse')
    CLim(gcm, [0 10])
    saveFigure(['/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH222_m2m_SC_inset',num2str(i)])
%     print(gcf, '-dpng', '-r300', ['/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH222_m2mean_SC_inset_',num2str(i),'.png']);
end
%% mis2mean-SG
[grainsSG_S,ebsd_SG_S.grainId,ebsd_SG_S.mis2mean] = calcGrains(ebsd_SG_S,'threshold',20*degree);
%%
close all;figure;
f=newMtexFigure('figsize','tiny');
[~,mP] =plot(ebsd_SG_S('ol'), ebsd_SG_S('ol').mis2mean.angle./degree,'micronbar','on');
mP.micronBar.length=100;
CLim(gcm, [0 10])
mtexColorbar('title', 'Misorientatation to mean orientation [º]')
hold on
plot(grainsSG_S('ol').boundary,'linewidth',1.5,'micronbar','on')
hold off
% 3 insets
xtd=ebsd_SG_S.extend;
xmid=round((xtd(2)-xtd(1))/7);
ymid=round((xtd(4)-xtd(3))/2);
sqSZ=100;

for i=[2,4,6]
    pos=[i*xmid-sqSZ ymid-sqSZ/2 sqSZ sqSZ];
    rectangle('Position',pos,'EdgeColor',[241 89 43]/255,'LineStyle','--','LineWidth',2)
end
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH222_mis2mean_SG')
mtexColorbar('off')
mP.micronBar.visible = 'off';
axis square tight
for i=[2,4,6]
    pos=[i*xmid-sqSZ ymid-sqSZ/2 sqSZ sqSZ];
    axis([pos(1) pos(1)+sqSZ pos(2) pos(2)+sqSZ])
%     az=180; el=-90; view(az, el);
    set(gca,'Ydir','reverse')
    CLim(gcm, [0 10])
    saveFigure(['/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH222_m2m_SG_inset',num2str(i)])
%     print(gcf, '-dpng', '-r300', ['/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/HH222_m2mean_SG_inset_',num2str(i),'.png']);
end
%% Pole figures at 7 GPa (HH221 and HH222)
hW=8*degree;
cs=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');

% HH221    

grains_SC_HH221=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH221/HH221-SC/map20190511174102326_crop_clean_grains.mat');
grains_SC_HH221=grains_SC_HH221.grains;

grains_SG_HH221=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH221/HH221-SG/map20190510212038753_crop_clean_grains.mat');   
grains_SG_HH221=grains_SG_HH221.grains;

odf_SC_HH221=calcDensity(grains_SC_HH221('ol').meanOrientation,'halfwidth',hW);
odf_SG_HH221=calcDensity(grains_SG_HH221('ol').meanOrientation,'halfwidth',hW);
% HH222    
grains_SC_HH222=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH222/HH222-SC/map20200825174846366_Rescan_clean_ps_grains.mat');
grains_SC_HH222=grains_SC_HH222.grains;

grains_SG_HH222=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/data/HH222/HH222-SG/map20200826214013815_Rescan_clean_ps_grains.mat');
grains_SG_HH222=grains_SG_HH222.grains;

odf_SC_HH222=calcDensity(grains_SC_HH222('ol').meanOrientation,'halfwidth',hW);
odf_SG_HH222=calcDensity(grains_SG_HH222('ol').meanOrientation,'halfwidth',hW);

%     
   h= Miller({1,0,0},{0,1,0},{0,0,1},cs,'hkl');
 
% pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'\sigma_3','\sigma_1'},...
%   'BackgroundColor','w','tag','axesLabels',varargin{:});
% setMTEXpref('pfAnnotations',pfAnnotations);
setMTEXpref('xAxisDirection','east');setMTEXpref('zAxisDirection','intoPlane'); %EDAX: A1 up(-y), A2 left(-x)
setMTEXpref('outerPlotSpacing',2); setMTEXpref('innerPlotSpacing',2);
setMTEXpref('figsize','medium')
figure
f=newMtexFigure('nRows',4,'ncols',3);

plotPDF(odf_SG_HH221, h,'lower','grid','on','grid_res',15*degree,'lineWidth',0.01)
f.nextAxis;
plotPDF(odf_SC_HH221, h,'lower','grid','on','grid_res',15*degree,'lineWidth',0.01)
f.nextAxis;
plotPDF(odf_SG_HH222, h,'lower','grid','on','grid_res',15*degree,'lineWidth',0.01)
f.nextAxis;
plotPDF(odf_SC_HH222, h,'lower','grid','on','grid_res',15*degree,'lineWidth',0.01)
mtexColorMap(parula)
CLim(gcm,'equal'); mtexColorbar ;% Single colorbar for all plots
f.drawNow
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/PF.png')
%% Pressure calibration
phaseTransition = {'','Bi_{I}-Bi_{II}','B_{II}-Bi_{III}','Bi_{III}-Bi_{IV}','ZnTe - ZnTe_{HP1}','ZnTe_{HP1} -ZnTe_{HP2}','ZnS-ZnS_m'};
pressLoad = [0, 0.472380952, 0.497777778, 1.320634921, 1.6, 2.85968254, 5.333333333]; %MN
samplePressure= [0, 2.55, 2.7, 7.7, 9.6, 12, 15.6];%GPa

fo = fitoptions('Method','NonlinearLeastSquares','Normalize','off','Algorithm','Trust-Region','Display','iter','TolFun',1e-15,'TolX',1e-15);

FitStr = 'a*exp(-x/b) + c';

[curve,~] = fit(pressLoad', samplePressure', FitStr, 'Robust', 'LAR'); % curve fitting
% Plot data and Fit
close all
figure('Units','normalized','pos',[0.1 0.1 0.45 0.35],'NumberTitle', 'off', 'Name', 'gif');

scatter(pressLoad, samplePressure,'filled')
hold on
LoadInterpValues=0:0.1:max(pressLoad);

plot(LoadInterpValues,[curve.a*exp(-LoadInterpValues/curve.b) + curve.c])
xlabel('Load [MN] '); ylabel('Pressure [GPa]');
axis square tight
box on
legend({'Pressure points', 'Fit'},'Location','northEastOutside')
% saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper4/illustrations/PressurePoints.pdf')