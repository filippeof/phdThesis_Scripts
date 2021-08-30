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

    How to cite:
    Ferreira, F., Hansen, L., & Marquardt, K. (2021). 
    The effect of grain boundaries on plastic deformation of olivine. Journal of Geophysical Research: Solid Earth, 126, e2020JB020273. 
    https://doi.org/10.1029/2020JB020273

%}
%% 
%obs:
% area345 (used in Fig. 1,4,7,8, 9 and 10)  corresponds to the files of section YZ stitched together
% ebsd3_5_6_7_8_11 (used in fig. 11) corresponds to the files of section XY stitched together
% Files are available at: https://doi.org/10.5281/zenodo.4793424
%
%E= Low-strain sections
%C= Moderate-strain sections
%G= High-strain sections
%% Figure 1:  IPF map

fname=  'AREA345_clean.ang';
[fPt,sampleName,~] = fileparts(fname);

ebsd =EBSD.load(fname,'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');

phaseName='olivine';
cs=ebsd(phaseName).CS;

minCI =0.01; % Minimum confidence index
minGS =20; % Minimum  grain size
minANG =2*degree; % Minimum  grain angle threshold
sF=5;% smooth factor for grains calculation (number of iterations)

close all; figure
oM = ipfHSVKey(cs);
oM.inversePoleFigureDirection = zvector;% 001 IPF colorcoding
% oM.colorPostRotation = reflection(yvector);% change legend to match TSL colors
color = oM.orientation2color(ebsd(phaseName).orientations);

setMTEXpref('figSize','large');
setMTEXpref('outerPlotSpacing',0);
plot(ebsdCleanRot,ebsdCleanRot.prop.iq)
mtexColorMap black2white
hold on; plot(ebsd(phaseName), color, 'FaceAlpha',0.9);% plot orientations
hold on; plot(grains.boundary,'k', 'linewidth', 0.01, 'edgealpha', 0.07);hold off;
az=0; el=90; view(az, el);
%% Figure 2: Microstructures 
    %Loop over the ebsd files and plot image quality:
    plot(ebsd,ebsd.iq)
%% Figure 3b:  Slip transfer schematics

figure
setMTEXpref('figSize','tiny');
a=linspace(0,pi/2,37);
angs=combvec(a,a)';
angle1=angs(:,1);
angle2=angs(:,2);

[X,Y] = meshgrid(angle1,angle2);
m_mtx=cos(X).*cos(Y);
h=pcolor(rad2deg(X),rad2deg(Y),m_mtx)
h.EdgeColor='none'
set(h, 'EdgeColor', 'none');
shading interp;
colorbar
axis('equal');axis('tight')
% title(' m'' Factor','fontSize',16)
mtexColorMap(parula)
h=mtexColorbar('fontSize',16);
ylabel(h, 'm'' Factor','fontSize',16)

xlabel('κ [ º ]','fontSize',16); ylabel('Ψ [ º ]','fontSize',16);
set(gca, 'fontSize',16)

set(gca,'XLim',[0 90])
set(gca,'XTick',(0:15:90))
set(gca,'YLim',[0 90])
set(gca,'YTick',(0:15:90))
caxis([0 1])
set(gca, 'Color', 'none')

% saveFigure('plots/mprime.png')

%% Figure 4:  Pole Figures
%pole figures of area345 file
h = Miller({1,0,0},{0,1,0},{0,0,1},cs)';
figure; plotPDF(ebsd_xy(phaseName).orientations,h)
%% Figure 5:  SPO plots
% Grain shape (SPO) PT0499 
ebsdE=load(E_YZ_fname); %Load section YZ, low strain
ebsdE=ebsdE.ebsd; %low strain 

ebsdC=load(C_YZ_fname);
ebsdC=ebsdC.ebsd;% moderate strain

ebsdG=load(G_YZ_fname);
ebsdG=ebsdG.ebsd;% high strain
% %Rotate EBSD data (small rotations to make the edges perfectly horizontal - crystal orientations are not important here)

ebsdE_rot=rotate(ebsdE, rotation('axis',zvector,'angle',175*degree));

ebsdC_rot=rotate(ebsdC, rotation('axis',zvector,'angle',180.5*degree));
% ebsdC_rot=rotate(ebsdC_rot, rotation('axis',yvector,'angle',180*degree),'center',mean([ebsdC_rot.prop.x,ebsdC_rot.prop.y]));

ebsdG_rot=rotate(ebsdG,rotation('axis',zvector,'angle',181*degree));
% ebsdG_rot=rotate(ebsdG_rot,rotation('axis',yvector,'angle',180*degree),'center',mean([ebsdG_rot.prop.x,ebsdG_rot.prop.y]));

% 
% figure
% plot(ebsdE_rot('ol'),ebsdE_rot('ol').orientations)
% az=0; el=90; view(az, el);
% 
% figure
% plot(ebsdC_rot('ol'),ebsdC_rot('ol').orientations)
% az=0; el=90; view(az, el);
% 
% figure
% plot(ebsdG_rot('ol'),ebsdG_rot('ol').orientations)
% az=0; el=90; view(az, el);
%% calc grains
sF=5;% smooth factor for grains calculation
hW = 8*degree; %ODF halfwidth
minCI =0.1; % Minimum confidence index
minGS =20; % Minimum  grain size
minANG =2*degree; % Minimum  grain angle threshold
misIndexing=false;
[ebsdE_rot,grainsE]=cleanGrains(ebsdE_rot,minCI,minANG,minGS,sF,misIndexing);
[ebsdC_rot,grainsC]=cleanGrains(ebsdC_rot,minCI,minANG,minGS,sF,misIndexing);
[ebsdG_rot,grainsG]=cleanGrains(ebsdG_rot,minCI,minANG,minGS,sF,misIndexing);
%
f=newMtexFigure;
plot(grainsE.boundary)
f.nextAxis;
plot(grainsC.boundary)
f.nextAxis;
plot(grainsG.boundary)
% plot rose diagram
ctE=mean([grainsE.x, grainsE.y])
ctC=mean([grainsC.x, grainsC.y])
ctG=mean([grainsG.x, grainsG.y])
%%
ttlE=['\gamma \approx ', num2str(1.8)];
ttlC=['\gamma \approx ', num2str(5.6)];
ttlG=['\gamma \approx ', num2str(9)];

sqSize=300;
sqHalfSize=sqSize/2;
grainsE_crop =grainsE(inpolygon(grainsE,[ctE(1)+sqHalfSize, ctE(2)+sqHalfSize, sqSize, sqSize]));
grainsC_crop =grainsC(inpolygon(grainsC,[ctC(1)+sqHalfSize, ctC(2)+sqHalfSize, sqSize, sqSize]));
grainsG_crop =grainsG(inpolygon(grainsG,[ctG(1)+sqHalfSize, ctG(2)+sqHalfSize, sqSize, sqSize]));

%
close all
grainsE_crop(grainsE_crop.grainSize<30)=[];
grainsC_crop(grainsC_crop.grainSize<30)=[];
grainsG_crop(grainsG_crop.grainSize<30)=[];

%
[omegaE, aE, bE] = grainsE_crop.fitEllipse;
[omegaC, aC, bC] = grainsC_crop.fitEllipse;
[omegaG, aG, bG] = grainsG_crop.fitEllipse;

%
setMTEXpref('figSize','medium');
f=newMtexFigure('nRows',1,'nCols',3);
plot(grainsE_crop,omegaE./degree,'lineColor','w')
mtexTitle(ttlE, 'Interpreter','tex')
mtexColorbar
% az=0; el=90; view(az, el);
CLim(gcm, [0 180])
mtexColorMap WhiteJet

%
f.nextAxis;
plot(grainsC_crop,omegaC./degree,'lineColor','w')
mtexTitle(ttlC, 'Interpreter','tex')
mtexColorbar
% az=0; el=90; view(az, el);
CLim(gcm, [0 180])
mtexColorMap WhiteJet

%
f.nextAxis;
plot(grainsG_crop,omegaG./degree,'lineColor','w')
mtexTitle(ttlG, 'Interpreter','tex')
mtexColorbar
% az=0; el=90; view(az, el);
CLim(gcm, [0 180])
mtexColorMap WhiteJet

%% Histograms
figure
subplot(1,3,1)
polarhistogram(omegaE, 25,'Normalization','probability');
thetalim([0 180]);
title(ttlE, 'Interpreter','tex')
ax=gca;
ax.RLim=[0 0.13];

subplot(1,3,2)
polarhistogram(omegaC, 25,'Normalization','probability');
thetalim([0 180]);
title(ttlC, 'Interpreter','tex')
ax=gca;
ax.RLim=[0 0.13];

subplot(1,3,3)
polarhistogram(omegaG, 25,'Normalization','probability');
thetalim([0 180]);
title(ttlG, 'Interpreter','tex')
ax=gca;
ax.RLim=[0 0.13];
%% Average Ellipses 
figure
plotEllipse([0,0],mean(aE),mean(bE),mean((pi/2)-omegaE),'lineColor','g','lineWidth',3,'DisplayName',[ttlE,', \theta = ',num2str(90-mean(omegaE./degree),'%.0f'), 'º'])
hold on
plotEllipse([0,0],mean(aC),mean(bC),mean((pi/2)-omegaC),'lineColor','y','lineWidth',3,'DisplayName',[ttlC,', \theta = ',num2str(90-mean(omegaC./degree),'%.0f'), 'º'])
hold on
plotEllipse([0,0],mean(aG),mean(bG),mean((pi/2)-omegaG),'lineColor','r','lineWidth',3,'DisplayName',[ttlG,', \theta = ',num2str(90-mean(omegaG./degree),'%.0f'), 'º'])
hold on
line([0, 0, nan, -10, 10]',[-10, 10, nan, 0, 0]','Color','k','HandleVisibility','off')
hold off
axis off
az=0; el=90; view(az, el);
hl=legend;
% set(hl, 'Interpreter','tex')

%% Figure 6:  GBPD
SegmentsFileList={...
 '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/E_xy.txt',
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/E_xz.txt',
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/E_yz.txt',
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/C_xy.txt',
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/C_xz.txt',
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/C_yz.txt',
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/G_xy2.txt',
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/G_xz.txt',
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/G_yz.txt'
};
%% Normalize segments from the minimum of every section (xy,xz,yz) and strain interval (E,C,G)
% Find the number of minimum segments of each file (function defined at the end of this file)
minLines = minNumberSegments(SegmentsFileList)
segProperty='numberOfSegments';

for ii=1:length(SegmentsFileList)
    segments_fname_in=SegmentsFileList{ii};
    [fPt,sampleName,~] = fileparts(segments_fname_in);
    segments_fname_out = ['/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/',sampleName,'_norm.txt'];
    filterSegments(segments_fname_in, segProperty,minLines,segments_fname_out);
end
%% join segments of starting material
SegmentsFileList_SM={...
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/map20181031224212864_segments.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/map20180811102351587_segments.txt'}
SM_segments_fname='/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/PT0535.txt';
joinSegments(SM_segments_fname, SegmentsFileList_SM)

%% Join segments of each section
SegmentsFileList_norm={...
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/E_xy_norm.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/E_xz_norm.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/E_yz_norm.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/C_xy_norm.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/C_xz_norm.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/C_yz_norm.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/G_xy2_norm.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/G_xz_norm.txt'
'/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/G_yz_norm.txt'
}

E_segmentsList=SegmentsFileList_norm(1:3)
E_segments_fname_norm='/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/E_norm.txt'
joinSegments(E_segments_fname_norm,E_segmentsList)

C_segmentsList=SegmentsFileList_norm(4:6)
C_segments_fname_norm='/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/C_norm.txt'
joinSegments(C_segments_fname_norm,C_segmentsList)

G_segmentsList=SegmentsFileList_norm(7:9)
G_segments_fname_norm='/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/G_norm.txt'
joinSegments(G_segments_fname_norm,G_segmentsList)
%% Separate segments in high and low angle gb segments
ECG_norm_list={...
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/E_norm.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/C_norm.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/G_norm.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/PT0535.txt'
    }

for ii=1:length(ECG_norm_list)
    segments_fname_in=ECG_norm_list{ii};

    [fPt,sampleName,~] = fileparts(segments_fname_in);
    segments_fname_out_lowAngle=['/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/',sampleName,'_LAGB.txt']; %** Change accordingly
    segments_fname_out_highAngle=['/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/',sampleName,'_HAGB.txt'];; %** Change accordingly

    segProperty='misAngle'; % Other properties: 'segLength','traceAngle','numberOfSegments'
    segValue_low=[0 20]; % low angle gb 
    segValue_high=[20 nan]; % high angle gb

    filterSegments(segments_fname_in, segProperty, segValue_low, segments_fname_out_lowAngle)
    filterSegments(segments_fname_in, segProperty, segValue_high, segments_fname_out_highAngle)
end
%% Plot GBPD
% plot list = all boundaries
                    %HAGB
                    %LAGB
gbpdPlotList={...
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/PT0535.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/E_norm.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/C_norm.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/G_norm.txt'
    
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/PT0535_HAGB.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/E_norm_HAGB.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/C_norm_HAGB.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/G_norm_HAGB.txt'
    
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/PT0535_LAGB.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/E_norm_LAGB.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/C_norm_LAGB.txt'
    '/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/data/segments_edax_2deg/normalizedSegments/G_norm_LAGB.txt'
    };

    cs = crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');
    GBCD_par=[5.0000, 1.0000, 1.0000, 9.0000, 9.0000, 0, 1.0000, 0]; %[msym, rot, rad, resolution, seg, stepSize, cutOff]
    GBCDfolder='/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/segmentsTest/GBPD-CD_MTEX-0.9.2/GBCD';
    setMTEXpref('figSize','small');
    figure; f=newMtexFigure('Name','GBPD', 'NumberTitle','off', 'nrows',3,'ncols',4);
    
for ii=1:length(gbpdPlotList)
    segFname=gbpdPlotList{ii};

    [fPt,sampleName,~] = fileparts(segFname);
    datFname = segments2dat(GBCDfolder, segFname, GBCD_par); 
    
    % Load vectors and m.u.d from the .dat file
    gbpd=loadGBPD(datFname); 

    % Plot GBPD
if ii~=1
    f.nextAxis;
end
    plotGBPD(gbpd,cs); 
    % title 
	numSegments =getNumberOfLines(segFname);
    ttl=mtexTitle(['n: ', num2str(numSegments)]);
    set(ttl, 'horizontalAlignment', 'right');
    set(ttl, 'units', 'normalized');
    ttlPos = get(ttl, 'position');
    set(ttl, 'position', [1.05 -0.1 ttlPos(3)])% bottom right
    
    %Plot Miller Index over GBPD
    h= Miller({1,0,0},{0,1,0},{0,0,1},cs,'hkl');
    hold on
    plot(h,'markerFaceColor',[0.5765 0.2314 0.2549],'markerSize',5, 'labeled','all', 'fontSize',14,'Color',[ 0.5765 0.2314 0.2549],'upper')
    hold off
    mtexColorbar('title','mud')
end
    %Add colorbar
    mtexColorMap(parula); 
%     CLim(gcm,'equal'); 
%     mtexColorbar('title','mud') ;% Single colorbar for all plots
%     print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/Illustrations/GBPD_rev_0421_max1-5mud.pdf','-dpdf', '-painters');
%% colorbar max 1.5
figure; f=newMtexFigure;
plot(vector3d.rand(100),'upper','contourf')
mtexColorbar
CLim(gcm, [0 1.5])
ax
% ax.delete;
print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/Illustrations/colorbar1-5mud.pdf','-dpdf', '-painters');


%% Figure 7:  Axis and angle of misorientation

%Axis/angle distribution
gb=grains.boundary(phaseName,phaseName);
[~,idxUnq]=unique(sortrows(gb.grainId,1),'rows'); % only one segment/ boundary
gb_mis=gb(idxUnq).misorientation;
% gb_low=gb(gb_mis.angle>3*degree & gb_mis.angle<20*degree);
% gb_high=gb(gb_mis.angle>20*degree);

% plot angle distribution
setMTEXpref('figSize','small');
    newMtexFigure('Name',sampleName,'NumberTitle','off');
    plotAngleDistribution(gb_mis, 'DisplayName', [phaseName,' - ', phaseName])
    hold on
    plotAngleDistribution(cs, cs, 'DisplayName', 'Uniform distribution')
    hold off
    saveFigure(fullfile(fPt,[sampleName,'_MAngD.png']))

% plot axis distribution
 setMTEXpref('figSize','small');
 mtexFig = newMtexFigure;
 mtexFig.nrows=3;

for ang=5:10:105
    gb_ang=gb(gb_mis.angle>=ang*degree & gb_mis.angle<(ang+10)*degree);
    plotAxisDistribution(gb_ang.misorientation, 'smooth', 'parent', mtexFig.nextAxis)
    mtexTitle([num2str(ang),' - ',num2str(ang+10),'^{\circ}'],'Interpreter','tex')
end
% mtexColorbar('multiple')
CLim(gcm,'equal'); % set equal color range to all plots
mtexColorbar 
saveFigure(fullfile(fPt,[sampleName,'_MAxisD.png']))

%% Figure 8:  Grain size and texture strength evolution as a function of strain
% plot mindex-grainSize evolution

pname = '/Users/fereira/Documents/DOC/GBCD_Olivine/Data/IRTGSeminar_seismicCoursePoster0917/YZ';
% which files to be imported
fname = [pname '/AREA345_clean.ang'];
% fname ='pt0499B_HR_clean.ang';
hold off; close all;%
addpath(genpath('Scripts')) % add functions folder to Path

CS =getCS(fname);
% create an EBSD variable containing the data
ebsd=EBSD.load(fname,'CS',CS,'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');
xtd=ebsd.extend;
ebsd=rotate(ebsd,rotation.byAxisAngle(zvector,180*degree),'center', [mean(xtd(1:2)), mean(xtd(3:4))]);
ebsd=rotate(ebsd,rotation.byAxisAngle(xvector,-180*degree),'center', [mean(xtd(1:2)), mean(xtd(3:4))]);
    
sF=2;% smooth factor for grains calculation
hW = 8*degree; %ODF halfwidth
minCI =0.1; % Minimum confidence index
minGS =20; % Minimum  grain size
minANG =2*degree; % Minimum  grain angle threshold

    cs=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine', 'color', 'green');
    phaseName=cs.mineral; % desired phase
%     ebsdCleanRot=rotate(ebsd,rotation('axis',zvector,'angle',(pi)));%,'KeepXY');
%     ebsdCleanRot=ebsd;
    [ebsdCleanRot,grains]=cleanGrains(ebsd,minCI,minANG,minGS,sF);
%
n=8; % number of ebsd sections
xtd=ebsdCleanRot.extend;
xmax=max(ebsdCleanRot.prop.x);
xstp=floor(xmax/n);
gamma=[];mGS=[];mindex=[];jindex=[];
for i=1:n
    pos=[(i-1)*xstp xtd(3) xstp xtd(4)];
    ind = inpolygon(ebsdCleanRot,pos); % select indices by rectangle
    ebsd_cut=ebsdCleanRot(ind);
%     [ebsd_cut,grains_cut]=cleanGrains(ebsd_cut,minCI,minANG,minGS,sF);
    grains_cut=grains(inpolygon(grains,pos));
    xmean=(i*xstp)-xstp/2;
    gamma(i)=(xmean/xmax)*10.9;
    mGS(i)=mean(grains_cut.diameter)
    odf = calcDensity(ebsd_cut(phaseName).orientations,'halfwidth',hW);
    mindex(i)=calcMIndex(odf);
    jindex(i)=textureindex(odf);
end
%%
setMTEXpref('figSize','medium');
colors=[ [0, 0.4470, 0.7410];...
            [0.8500, 0.3250, 0.0980];...
            [0.9290, 0.6940, 0.1250];...
            [0.4940, 0.1840, 0.5560];...
            [0.4660, 0.6740, 0.1880];...
            [0.6350, 0.0780, 0.1840];...
            [0.3010, 0.7450, 0.9330] ];
        
newMtexFigure;
subplot(1,2,1)
yyaxis left
xx = min(gamma):.1:max(gamma);
yy = spline(gamma,mGS,xx);
plot(xx,yy,'k','LineWidth', 4,'Color',[colors(1,:),0.5],'DisplayName', 'Average grain size [µm]')

hold on 
% scatter(gamma,mGS,10,colors(1,:),'filled','HandleVisibility','off')
err=ones(size(mGS));
errorbar(gamma,mGS,err,'.','Color',[colors(1,:),0.5],'HandleVisibility','off')

hold off
ylabel('Grain Size [µm]');
xlim([0 10.9])
axis square 
set(gca,{'ycolor'},{'k'})  % Left color red, right color blue...

yyaxis right
xx = min(gamma):.1:max(gamma);
yy = spline(gamma,mindex,xx);
plot(xx,yy,'LineWidth', 4,'Color',[colors(2,:),0.5],'DisplayName', 'M-index')
hold on 
scatter(gamma,mindex,10,colors(2,:),'filled','HandleVisibility','off')
hold off
ylabel('M-index');
xlabel('Strain [\gamma]');
hold off   
xlim([0 10.9])
axis square 

set(gca,{'ycolor'},{'k'})  % Left color red, right color blue...
legend

subplot(1,2,2)
yyaxis right
xx = min(gamma):.1:max(gamma);
yy = spline(gamma,jindex,xx);
plot(xx,yy,'LineWidth', 4, 'Color',[colors(3,:),0.5], 'DisplayName', 'J-index')
hold on 
scatter(gamma,jindex,10,colors(3,:),'filled','HandleVisibility','off')
hold off
ylabel('J-index');
axis square 
xlim([0 10.9])
axis square 
set(gca,{'ycolor'},{'k'})  % Left color red, right color blue...
legend


%% Figure 9:  m' factor map
% plot gb colored with m' factor

[ebsd,grains]=cleanGrains(ebsd,minCI,minANG,minGS,sF);
[gb,mp,rbv,grains]= slipTrev(grains,phaseName,'xy');
setMTEXpref('outerPlotSpacing',4); setMTEXpref('innerPlotSpacing',0);

        figure; plot(gb,mp,'lineWidth',0.4)
        az=180; el=-90; view(az, el);
        caxis([0 1])
        mtexColorbar ('title','{\it m''} factor')
        print(gcf, '-dtiff', '-r450', ['/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/Illustrations/',sampleName,'_mprime_rev.tiff']);
% make insets
setMTEXpref('figSize','small');
        figure; plot(gb,mp,'lineWidth',1.5)

xtd=ebsd.extend;
xmax=max(ebsd.prop.x);
xmid=floor((xtd(2)/7)/2);
ymid=xtd(4)/2;

for i=[2,7,12]
    pos=[i*xmid-125 ymid-125 250 250];
%     rectangle('Position',pos,'LineStyle','--','EdgeColor','k','lineWidth',2)
    rectangle('Position',pos,'LineStyle','-','EdgeColor','k','lineWidth',1)

    axis([pos(1) pos(1)+251 pos(2) pos(2)+251])
%     set(gca,'Ydir','reverse')
    az=180; el=-90; view(az, el);
    
    print(gcf, '-dtiff', '-r300', ['/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/Illustrations/',sampleName,'_inset',num2str(i),'_rev.tiff']);
end
    axis equal tight
% az=0; el=90; view(az, el);
        mtexColorbar ('title','{\it m''} factor')

  saveFigure(['/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/Illustrations/',sampleName,'insetPositions_rev.png'])

%% Figure 10:  Slip system activity

% Define Slip systems [n; b]
slipS(:,:,1)=[ 0 1 0  ; 1 0 0]; % a type
slipS(:,:,2)=[ 0 1 0  ; 0 0 1]; % b-type
slipS(:,:,3)=[ 1 0 0  ; 0 0 1]; % c-type
slipS(:,:,4)=[ 0 1 1  ; 1 0 0]; % d-type
slipS(:,:,5)=[ 0 0 1  ; 1 0 0]; % e type
slipS(:,:,6)=[ 1 1 0  ; 0 0 1];% add new ss
slipS(:,:,7)=[ 1 3 0  ; 0 0 1];% add new ss
% slipS(:,:,8)=[ n1 n2 n3  ; b1 b2 b3];% add new ss

nSS=size(slipS,3); % number of slip systems

cs=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');
CRSS=ones(nSS,1);% Same CRSS
clear n b SS
for i =1:nSS
    n{i}=Miller(slipS(1,1,i),slipS(1,2,i),slipS(1,3,i),cs,'hkl');% n=slip plane normal
    b{i}=Miller(slipS(2,1,i),slipS(2,2,i),slipS(2,3,i),cs,'uvw');%b=slip direction
    % check if angle b to n is 90 degrees (i.e. that slip direction is in the slip plane)
    Angle_n_to_b = round(angle(n{i},b{i})./degree);
    if Angle_n_to_b ==90
        SS{i} = slipSystem(b{i},n{i},CRSS(i));
    end
end
SS=[SS{:}];
%
colors=[ [0, 0.4470, 0.7410];...
            [0.8500, 0.3250, 0.0980];...
            [0.9290, 0.6940, 0.1250];...
            [0.4940, 0.1840, 0.5560];...
            [0.4660, 0.6740, 0.1880];...
            [0.6350, 0.0780, 0.1840];...
            [0.3010, 0.7450, 0.9330] ];
        %if more than 7 ss add new color 
%
legend_str={};
close all
figure('Position',[0 0 380 280])
g_ol=grains('olivine');
n=8; %number of points to interpolate
g_ctd=g_ol.centroid;
maxX=max(g_ctd(:,1));
gamma=( (maxX-g_ctd(:,1)) /maxX)*10.9;
step_gamma=10.9/n;
proport_slip_i=zeros(nSS,n);
middle_gamma=zeros(1,n);
activessId=g_ol.prop.activeSSid_sym;
for j=1:n
    id=(gamma>(j-1)*step_gamma) & (gamma<=j*step_gamma);
    middle_gamma(j)=mean([(j-1)*step_gamma,j*step_gamma]);

    for i=1:nSS
        proport_slip_i(i,j)= abs(length(find(activessId(id)==i))/length(activessId(id)));
    end
end

for i=1:nSS
    xx = min(middle_gamma):.1:max(middle_gamma);
    yy = spline(middle_gamma,proport_slip_i(i,:),xx);
    plot(xx,yy,'-','LineWidth',4,'Color',[colors(i,:),0.5]);
    % plot (middle_gamma, proport_slip_i(i,:))
    legend_str{i}= [' [',num2str(round(SS(i).b.uvw)),']', ' (',num2str(round(SS(i).n.hkl)) ,')' ];
    hold on
    scatter(middle_gamma,proport_slip_i(i,:),10,colors(i,:),'filled','HandleVisibility','off')
end
% hold off
 legend_str = legend_str(~cellfun('isempty',legend_str)); % clear empty legend
        
% scatter gamma x m'
maxX=max(gb.midPoint(:,1));
gamma=( (maxX-gb.midPoint(:,1)) /maxX)*10.9;
step_gamma=10.9/n;
mean_mp=[]; std_mp=[]; SE=[]; middle_gamma=[];
id_del=(gb.midPoint(:,1)>1940 & gb.midPoint(:,1)<1960);
gb(id_del)=[];
mp(id_del)=[];
gamma(id_del)=[];

for i=1:n
    id=(gamma>(i-1)*step_gamma) & (gamma<=i*step_gamma);
    mean_mp(i)=mean(mp(id));
    std_mp(i)=std(mp(id));
    SE(i)= std(mp(id))/sqrt(length(mp(id)));
    middle_gamma(i)=mean([(i-1)*step_gamma,i*step_gamma]);
end

% figure
% errorbar(middle_gamma,mean_mp,std_mp) 
% hold on
% scatter(middle_gamma,mean_mp)
% hold off

hold on
xx = min(middle_gamma):.1:max(middle_gamma);
yy = spline(middle_gamma,mean_mp,xx);
scatter(middle_gamma,mean_mp,10,[0 0 0],'filled','HandleVisibility','off')
hold on
plot(xx,yy,'Color',[0,0,0,0.5],'LineWidth', 4)
hold off

legend_str{end+1}='{\it m''} factor';
l=legend(legend_str, 'fontSize',10,'Location','northEastOutside');

axis square tight
axis([0 10.9 -0.001 1])
xlabel('Strain [{\it\gamma} ]'); ylabel('Inferred slip-system activity');
set(gca,'fontSize',12)
g=gca; 
% 
box on
ax2= axes('Position',g.Position, 'color', 'none',...
'XTick',[], 'YAxisLocation', 'right')
ylabel(ax2, '{\it m''} factor');
axis square tight
axis([0 10.9 0 1])
 set(gca,'fontSize',12)
print(['/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/Illustrations/',sampleName,'_mprime_ssXstrain'],'-dpdf', '-painters') 


%% Figure 11: KAM grid 

ebsd.id=(1:length(ebsd))';

grains=smooth(grains,5);

%
gb_all=grains.boundary('ol','ol');
gb_all=[gb_all,grains.innerBoundary];
gb_step=1;% Sample every 'gb_step' grain boundary
gb=gb_all(1:gb_step:end);

%
gbdir=gb.direction;
trace_dir=angle(gbdir,xvector);

% boundary mid-point (x,y)
xy_mid = gb.midPoint;
% 
% % xy of gb traces

smpSZ=30; % number of pixels from gb (smpSZ/2 for each side)

X1_o = xy_mid(:,1);
Y1_o = xy_mid(:,2);

X2_o = xy_mid(:,1) - 0.5 * smpSZ .* sin(trace_dir);
Y2_o = xy_mid(:,2) + 0.5 * smpSZ .* cos(trace_dir);

X3_o = xy_mid(:,1) + 0.5 * smpSZ .* sin(trace_dir);
Y3_o = xy_mid(:,2) - 0.5 * smpSZ .* cos(trace_dir);


ind=((gbdir.x>0 & gbdir.y<0) | (gbdir.x<0 & gbdir.y>0));

X2_o(ind) = xy_mid(ind,1) - 0.5 * smpSZ .* sin(pi-trace_dir(ind));
Y2_o(ind) = xy_mid(ind,2) + 0.5 * smpSZ .* cos(pi-trace_dir(ind));

X3_o(ind) = xy_mid(ind,1) + 0.5 * smpSZ .* sin(pi-trace_dir(ind));
Y3_o(ind) = xy_mid(ind,2) - 0.5 * smpSZ .* cos(pi-trace_dir(ind));

btx = reshape([X1_o,X2_o,nan(size(X1_o))].',[],1);
bty = reshape([Y1_o,Y2_o,nan(size(Y1_o))].',[],1);
    
btx2 = reshape([X1_o,X3_o,nan(size(X1_o))].',[],1);
bty2 = reshape([Y1_o,Y3_o,nan(size(Y1_o))].',[],1);

X1_o=[X1_o; X1_o]; %midpoint x
Y1_o=[Y1_o; Y1_o];%midpoint y
X2_o=[X2_o; X3_o]; %end point x
Y2_o=[Y2_o; Y3_o];%end point y
    
%   Plot profiles from/to GB:

% figure
% kam_all=KAM(ebsd,'threshold',10*degree);
% plot(ebsd,kam_all/degree,'figSize','huge')
% hold on
% 
% plot(gb_all,'linewidth',0.5,'color','k')
% hold on
% plot(btx,bty,'-','linewidth',1.5,'color','m','DisplayName','boundary traces');
% hold on
% plot(btx2,bty2,'-','linewidth',1.5,'color','c','DisplayName','boundary traces');
%  
%     hold off
%     axis equal tight
%     mtexColorbar('title','KAM [�]');
%
%     legend('show')
%     zoom(30)
%

% GET kam of [[X1_o Y1_o];[X2_o Y2_o]]


kam_all=KAM(ebsd,'threshold',10*degree,'order',3);
kam=[];

mis2mean_all=ebsd.mis2mean.angle;
mis2mean_ol=[];

dist=[]; 

f = waitbar(0,'1','Name','Calculating... ',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

for i=1:1000:length(X1_o)
    lineSec =  [X1_o(i) Y1_o(i); X2_o(i) Y2_o(i)];
    [ebsd_line, d] = spatialProfile(ebsd,lineSec);
    kam=[kam; kam_all(ebsd_line.id)];
    mis2mean_ol=[mis2mean_ol; mis2mean_all(ebsd_line.id)];
    dist=[dist; d];
    waitbar(i/length(X1_o),f,sprintf('%d',i))
     if getappdata(f,'canceling')
        break
    end
end
delete(f)
mis2mean_ol=mis2mean_ol/degree;
kam=kam/degree;
%Scatter plot

% 
% figure
% scatter(dist,kam)
% xlabel('Distance from gb (um)'); ylabel('KAM [�]');
% 
% figure
% scatter(dist,mis2mean_ol)
% xlabel('Distance from gb (um)'); ylabel('Mis2mean [�]');
%% 2D density plot
close all
setMTEXpref('figSize','tiny');
ind_kam=kam>1.5 & kam<8;

dat_kam= [dist(ind_kam),kam(ind_kam)];
% plot kam
f=newMtexFigure;
n=hist3(dat_kam,[smpSZ smpSZ]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat_kam(:,1)),max(dat_kam(:,1)),size(n,1)+1);
yb = linspace(min(dat_kam(:,2)),max(dat_kam(:,2)),size(n,1)+1);

pcolor(xb,yb,n1,'parent',f.gca);
xlabel('Distance from gb [um]'); ylabel('KAM [�]');
mtexColorbar('title','Frequency') ;
% shading interp
% set(gca,'ColorScale','log');
axis square tight
saveFigure('/Users/fereira/Documents/MATLAB/pt04993d/gb_proximitkambiggerthan2.png')
print( '/Users/fereira/Documents/MATLAB/pt04993d/gb_proximitkambiggerthan2','-dpdf', '-painters')
%% plot mis 2 mean
ind_mis2mean=mis2mean_ol>1.5 & mis2mean_ol<20 ;
dat_mis2mean= [dist(ind_mis2mean),mis2mean_ol(ind_mis2mean)];

f2=newMtexFigure;
n=hist3(dat_mis2mean,[10 10]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;

xb = linspace(min(dat_mis2mean(:,1)),max(dat_mis2mean(:,1)),size(n,1)+1);
yb = linspace(min(dat_mis2mean(:,2)),max(dat_mis2mean(:,2)),size(n,1)+1);

pcolor(xb,yb,n1,'parent',f2.gca);
xlabel('Distance from gb [um]'); ylabel('mis2mean [�]');
mtexColorbar('title','Frequency') ;
% shading interp

axis square tight
saveFigure('/Users/fereira/Documents/MATLAB/pt04993d/gb_proximitmis2mean_olbiggerthan2.png')
print( '/Users/fereira/Documents/MATLAB/pt04993d/gb_proximitmis2mean_olbiggerthan2','-dpdf', '-painters')
%% profile highliting GB
fname='/Users/fereira/Documents/MATLAB/pt04993d2/map20180805215930583.ang';

ebsd = loadEBSD(fname,'interface','ang',...
  'convertEuler2SpatialReferenceFrame');
ebsd('ol').CS=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');
       
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'threshold',20*degree);
ebsd(grains(grains.grainSize<50)) = [];
ebsd('En')=[];
ebsd.id=1:length(ebsd);
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'threshold',20*degree);
%
grains=smooth(grains,5);

% n=length(grains);
% randcolors=[rand(n,1),rand(n,1),rand(n,1)];
% plot(grains,randcolors,'figSize','huge')
kam_all=KAM(ebsd,'threshold',10*degree);
plot(ebsd,kam_all,'figSize','huge')
hold on;
plot(grains.boundary,'lineColor','w','lineWidth',2)
h = imline;
line_pos = wait(h);
delete(h)
hold on
line(line_pos(:,1),line_pos(:,2),'Color','m','LineWidth',4)
[ebsd_line,dist,idList]  = spatialProfile(ebsd,line_pos);
kam_line=kam_all(ebsd_line.id);

[g_id_crossed,~,~] = unique(ebsd_line.grainId,'stable');
% g_selected=grains(grains.id==g_id_crossed(1));
% for i =2:length(g_id_crossed)
%     g_selected=[g_selected,grains(grains.id==g_id_crossed(i))];
% end
g_selected=grains(g_id_crossed);
%  hold on
% scatter(ebsd_line.prop.x,ebsd_line.prop.y)
% hold off

hold on;
plot(g_selected.boundary,'lineColor','c','lineWidth',2)
hold off

g_ctd=g_selected.centroid;
text(g_ctd(:,1),g_ctd(:,2),num2str(g_selected.id))
%


figure('Units','normalized','pos',[0 0 1 0.5]);
changedIndexes = diff(ebsd_line.grainId)~=0;

gb_positions=dist(changedIndexes);
gb_positions=[gb_positions';gb_positions';nan(1,length(gb_positions))];
gb_positions=gb_positions(:);

% scatter(dist,ebsd_line.mis2mean.angle./degree)
% scatter(dist,kam_line/degree,'DisplayName', 'Pixel KAM');
yyaxis left
plot(dist,kam_line/degree,'DisplayName', 'Pixel KAM');
ylabel('KAM [�]');

yyaxis right
plot(dist,ebsd_line.mis2mean.angle./degree,'DisplayName', 'Pixel mis2mean');
ylabel('Misorientation to mean[�]');

y_gb=[ylim,nan]';
y_gb=repmat(y_gb,length(gb_positions)/3,1);
hold on
plot(gb_positions,y_gb,'--k','DisplayName', 'Grain boundaries')
hold off
legend('Location','northEastOutside')
axis tight
xlabel('Distance [um]'); 

%% Figure 13:  Grain boundary formation - Stereographic projections
%GB formation as result of disCreep: Tilt boundaries
% this code works with mtex 5.2.1 but did not work with 5.4
close all; clear all
% Define Slip systems [n; b]
slipS(:,:,1)=[ 0 1 0  ; 1 0 0];% a type
slipS(:,:,2)=[ 0 1 0  ; 0 0 1];% b-type
slipS(:,:,3)=[ 1 0 0  ; 0 0 1];% c-type
slipS(:,:,4)=[ 0 1 1  ; 1 0 0];% d-type
slipS(:,:,5)=[ 0 0 1  ; 1 0 0];% e type
slipS(:,:,6)=[ 1 1 0  ; 0 0 1];% add new ss
slipS(:,:,7)=[ 1 3 0  ; 0 0 1];% add new ss
% slipS(:,:,6)=[ n1 n2 n3  ; b1 b2 b3];% add new ss

nSS=size(slipS,3); % number of slip systems

cs=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');
CRSS=ones(nSS,1);% Same CRSS
clear n b SS
for i =1:nSS
    n{i}=Miller(slipS(1,1,i),slipS(1,2,i),slipS(1,3,i),cs,'hkl');% n=slip plane normal
    b{i}=Miller(slipS(2,1,i),slipS(2,2,i),slipS(2,3,i),cs,'uvw');%b=slip direction
    % check if angle b to n is 90 degrees (i.e. that slip direction is in the slip plane)
    Angle_n_to_b = round(angle(n{i},b{i})./degree);
    if Angle_n_to_b ==90
        SS{i} = slipSystem(b{i},n{i},CRSS(i));
    end
end
SS=SS';
n=[n{:}]; b=[b{:}];SS=vertcat(SS{:});
SS_sym= SS.symmetrise;
dS = dislocationSystem(SS);

%
setMTEXpref('outerPlotSpacing',40);
setMTEXpref('innerPlotSpacing',10);
% setMTEXpref('xAxisDirection','east');setMTEXpref('zAxisDirection','outOfPlane');
setMTEXpref('bAxisDirection','north');
setMTEXpref('aAxisDirection','east');


figure
f=newMtexFigure('nrows',3, 'ncols',3);
h = Miller({1,0,0},{0,1,0},{0,0,1},cs)';
ang=[5:2:85]*degree;

for i=1:length(SS)
%     b_nice=[SS(i).b.hkl]./norm([SS(i).b.hkl]); % slip direction
%     n_nice=[SS(i).n.hkl]./norm([SS(i).n.hkl]); % slip plane
    b_nice=SS(i).b.uvw; % slip direction
    n_nice=SS(i).n.hkl; % slip plane

    txt=sprintf('[%u%u%u] (%u%u%u)',round([b_nice, n_nice]));

    
    rot = rotation.byAxisAngle(dS(i).l, ang); %rotation axis = line vector
    h_rot=symmetrise(round(Miller(rot*(SS(i).b),cs)));
    maxIdx=5;

    m=round(h_rot); m=unique(m);
    m=symmetrise(m,cs,'upper');m=m(:);
    m=m(m.h<=maxIdx & m.k<=maxIdx & m.l<=maxIdx & sum(abs(m.hkl),2)<(2*maxIdx));


    if i~=1
        f.nextAxis;
    end

    plot(m,'upper','grid','grid_res',5*degree,'all','labeled','projection','eangle','fundamentalRegion','FontSize',9,'Interpreter','tex')
    mtexTitle(txt,'Interpreter','tex','FontSize',12)

end
f.drawNow;
% saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/Illustrations/GB_formation_ss_eangle.png')
print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper2/Illustrations/GB_formation_ss_eangle','-dpdf', '-painters')

%%
function minLines = minNumberSegments(segFnameList)
% Finds the minimum number of lines among all files in the list (segFnameList)
% Header lines (commented with #) are not counted
        % Normalize by minimum of segments among all files
        numlines=zeros(1,length(segFnameList));
        % Reads the number of lines of each file
       for j=1:length(segFnameList)
            A = regexp(fileread(segFnameList{j}), '\n', 'split')';
            A=cellfun(@strtrim, A,'UniformOutput' ,false); %remove empty lines
            A=A(~cellfun('isempty',A)); %remove empty lines

            %Number of header lines
            hL=0;
            for nline=1:50 %check first 50 lines
                if ischar(A{nline}) && startsWith(A{nline},'#')
                  hL = hL+1;
                end
            end
            numlines(j)=length(A)-hL;
        end
        %
        %file with minimum number of lines
        minLines=min(numlines);
end
%%
function numlines =getNumberOfLines(fileName)
    A = regexp(fileread(fileName), '\n', 'split')';
    A=cellfun(@strtrim, A,'UniformOutput' ,false); %remove empty lines
    A=A(~cellfun('isempty',A)); %remove empty lines

    %Number of header lines
    hL=0;
    for nline=1:50 %check first 50 lines
        if ischar(A{nline}) && startsWith(A{nline},'#')
          hL = hL+1;
        end
    end
    numlines=length(A)-hL;
end

function [ebsd,grains]=cleanGrains(ebsd,minPatternQuality,minANG,minGS,sF)
    % clean grains grains
    if isfield(ebsd.prop,'ci')
        ebsd(ebsd.prop.ci<minPatternQuality).phase=-1; 
    elseif isfield(ebsd.prop,'bc')
        ebsd(ebsd.prop.bc<minPatternQuality).phase=-1; 
    else
        warning('No EBSD quality filter')
    end

    [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'threshold',minANG);

    %Remove holes
    try
        notIndexed = grains('notIndexed');
        toRemove = notIndexed(notIndexed.grainSize ./ notIndexed.boundarySize<0.8);
        ebsd(toRemove) = [];
    catch
    end
    % done in OIM:
    % %Correct misindexing (only for shape analyses.. what to do about
    % orientation?)
    % for i=1:size(misIndexing.phase,1)
    %     gB = grains.boundary(misIndexing.phase,misIndexing.phase);
    %     rot = rotation('axis',misIndexing.axis,'angle',misIndexing.angle);
    %     ind = angle(gB.misorientation,rot)<misIndexing.tolerance;
    %     grains= merge(grains,gB(ind));
    % end
        ebsd(grains(grains.grainSize<minGS)) = [];

    % [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'threshold',minANG,'boundary','tight');
    [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'threshold',minANG);
    if sF~=0
        grains=smooth(grains,sF);
    end
end