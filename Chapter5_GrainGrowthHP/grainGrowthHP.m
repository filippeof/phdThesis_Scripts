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
%%
% Grain growth at HP - Grain growth Histograms

% 1 - Save grains data as mat file for all files

 % Parameters
 basePath_6 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/6pctEn';
basePath_13 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/13pctEn';
 allPaths={
    fullfile(basePath_6,'startingMaterial')       % Starting material FSG4 (FSG4-GMF2)
    fullfile(basePath_13,'startingMaterial')	  % Starting material FSG5
    fullfile(basePath_6,'10gpa/0C-0h')           % Z1993 (MA3B) - No heating

    fullfile(basePath_6,'5gpa/1400C-24h')      % Z1962
    fullfile(basePath_6,'7gpa/1400C-24h')      % Z1965
    fullfile(basePath_6,'12gpa/1400C-24h')    % Z1968
    
    fullfile(basePath_13,'1gpa/1400C-24h')    % B1272
    fullfile(basePath_13,'1gpa/1200C-72h')    % B1273
    fullfile(basePath_13,'1gpa/1200C-24h')    % B1274
    
    fullfile(basePath_13,'1gpa/1050C-24h')    % B1275
    fullfile(basePath_13,'1gpa/1200C-12h')    % B1276
    fullfile(basePath_13,'10gpa/1400C-24h')  % Z2032
    
    fullfile(basePath_13,'7gpa/1200C-24h')    % Z2033
    fullfile(basePath_13,'7gpa/1200C-12h')    % Z2034
    fullfile(basePath_13,'7gpa/1050C-24h')    % Z2035
    
    fullfile(basePath_13,'7gpa/1400C-12h')    % Z2047
    fullfile(basePath_13,'7gpa/1520C-24h')    % Z2049
    fullfile(basePath_13,'7gpa/1400C-72h')    % Z2051
    
    fullfile(basePath_13,'7gpa/1400C-24h')    % Z2060
    fullfile(basePath_13,'1gpa/1400C-12h')    % A1178
    fullfile(basePath_13,'1gpa/1400C-72h')    % A1179
    
    fullfile(basePath_13,'7gpa/1400C-48h')    % Z2062
    fullfile(basePath_13,'1gpa/1400C-8h')      % A1182
    
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-StartingMaterial')    % HV806-Starting Material
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-HighPressure')      % HV806- High Pressure
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-HighTemperature')    % HV806- High Temperature
    };
 for nbP=1:length(allPaths)
     basePath=allPaths{nbP};
    ebsdFileList= listEBSD(basePath,'.ang');
    fEN=[];

    % Cleaning parameters
    % Minimum confidence index
    minCI =0.01; 
    % Minimum  grain size (points)
    minGS =20; 
    % Minimum grain misorientation angle threshold
    minANG =20*degree; 
    % Grain boundaries smooth factor (number of iterations)
    sF=5;
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

    % Loop over all ebsd data, save grains data and plot maps      
    for i=1:length(ebsdFileList)
        ebsdFname=ebsdFileList{i};
        [fPt,sampleName,~] = fileparts(ebsdFname);
        grains_fname=fullfile(fPt,[sampleName,'_grains.mat']);
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
         % Smooth grain boundaries
        grains=smooth(grains,sF);
        % Remove outer boundary grains
        outerBoundary_id = any(grains.boundary.grainId==0,2);
        grain_id = grains.boundary(outerBoundary_id).grainId;
        grain_id(grain_id==0) = [];
        grains('id', grain_id) = [];

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
        
        hold on; plot(grains, 'FaceAlpha',0.85);% plot orientations
        hold on; plot(grains('indexed').boundary,'k', 'linewidth', 0.9, 'edgealpha', 0.6);hold off;
        lgd=legend; lgd.String{1}=['Olivine (',num2str(round( sum(grains('ol').area) /sum(grains.area) *100)),' %)'];
        lgd.String{2}=['Enstatite (',num2str(100-round( sum(grains('ol').area) /sum(grains.area) *100)),' %)'];
        set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[.9;.9;.9;.75]));  % [.9,.9,.9] is light gray; 0.75 means 25% transparent
        saveFigure(fullfile(fPt, [sampleName,'_phaseMap.png']))

        % ------- Save data-----------
        save (grains_fname, 'grains');
        fEN=[fEN,1-(sum(grains('ol').area) /sum(grains.area))];
    end
 end
%% Save grain size data for all samples
phaseName='olivine';
method='grid';
showPlot=false;
basePath_6 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/6pctEn';
basePath_13 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/13pctEn';

headers={'P(GPa)','T(ºC)','t(hours)','meanDiameter','meanInterceptLength','modeFit','muFit','sigmaFit','EnstatiteFraction', 'NumberOfGrains'};
filename=fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/AlldataSummary.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'%s \t', headers{:});
fprintf(fileID,'\n');
fclose(fileID);

exp_fname={

    fullfile(basePath_6,'startingMaterial')       % Starting material FSG4 (FSG4-GMF2)
    fullfile(basePath_13,'startingMaterial')	  % Starting material FSG5
    fullfile(basePath_6,'10gpa/0C-0h')           % Z1993 (MA3B) - No heating

    fullfile(basePath_6,'5gpa/1400C-24h')      % Z1962
    fullfile(basePath_6,'7gpa/1400C-24h')      % Z1965
    fullfile(basePath_6,'12gpa/1400C-24h')    % Z1968
    
    fullfile(basePath_13,'1gpa/1400C-24h')    % B1272
    fullfile(basePath_13,'1gpa/1200C-72h')    % B1273
    fullfile(basePath_13,'1gpa/1200C-24h')    % B1274
    
    fullfile(basePath_13,'1gpa/1050C-24h')    % B1275
    fullfile(basePath_13,'1gpa/1200C-12h')    % B1276
    fullfile(basePath_13,'10gpa/1400C-24h')  % Z2032
    
    fullfile(basePath_13,'7gpa/1200C-24h')    % Z2033
    fullfile(basePath_13,'7gpa/1200C-12h')    % Z2034
    fullfile(basePath_13,'7gpa/1050C-24h')    % Z2035
    
    fullfile(basePath_13,'7gpa/1400C-12h')    % Z2047
    fullfile(basePath_13,'7gpa/1520C-24h')    % Z2049
    fullfile(basePath_13,'7gpa/1400C-72h')    % Z2051
    
    fullfile(basePath_13,'7gpa/1400C-24h')    % Z2060
    fullfile(basePath_13,'1gpa/1400C-12h')    % A1178
    fullfile(basePath_13,'1gpa/1400C-72h')    % A1179
    
    fullfile(basePath_13,'7gpa/1400C-48h')    % Z2062
    fullfile(basePath_13,'1gpa/1400C-8h')      % A1182
    
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-StartingMaterial')    % HV806-Starting Material
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-HighPressure')      % HV806- High Pressure
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-HighTemperature')    % HV806- High Temperature
  
};

sampleNames = {'FSG4-SM','FSG5-SM','Z1993','Z1962','Z1965','Z1968','B1272','B1273', 'B1274'...
    ,'B1275','B1276', 'Z2032', 'Z2033', 'Z2034', 'Z2035','Z2047', 'Z2049',...
    'Z2051','Z2060', 'A1178', 'A1179', 'Z2062','A1182', 'HV806-SM', 'HV806-HP', 'HV806-HT'};
P_exp=[0.7, 0.7, 10, 5, 7, 12, 1, 1, 1, 1, 1, 10, 7, 7, 7, 7, 7, 7 7, 1, 1, 7, 1, 1E-11, 1, 1E-4];  % Pressure in GPa
T_exp=[1200, 1200, 25, 1400, 1400, 1400, 1400, 1200, 1200, 1050, 1200, 1400, 1200, 1200, 1050, 1400, 1520, 1400, 1400, 1400, 1400, 1400, 1400, 1250, 1400, 1400]; % Temperature in Cº
t_exp=[2, 2, 0.5, 24, 24, 24, 24, 72, 24, 24, 12, 24, 24, 12, 24, 12, 24, 72, 24, 12, 72, 48, 8, 3, 24, 24]; % Time in hours

data=nan(numel(headers),length(exp_fname));

for ii=1:length(exp_fname)
    tic
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat')
    nGrains=0; areaEn=[]; areaFo=[];
    phDiameter=[]; accumulatedSegLength=[];
    parfor jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        try
            areaEn = [areaEn, sum(grains('En').area)];
        catch
             areaEn = [areaEn, 0];
        end
        areaFo = [areaFo, sum(grains(phaseName).area)];
        nGrains=nGrains+length(grains(phaseName));
        phDiameter=[phDiameter; grains(phaseName).diameter];
        nLines=round(0.05*length(grains));
        if nLines<20; nLines=20; end
        accumulatedSegLength=[accumulatedSegLength, linearInterceptSP(grains,phaseName,method,nLines,showPlot)];  
    end
    fEn=sum(areaEn)/(sum(areaEn)+sum(areaFo)) ;
    meanDiameter=mean(phDiameter);
    meanInterceptLength=nanmean(accumulatedSegLength);
    pd= fitdist(phDiameter,'lognormal'); 
    mu=pd.mu;sigma=pd.sigma;
%     m_d=exp(mu+0.5*sigma^2);% mean
%     md_d=exp(mu); % median
    mod_d=exp(mu-sigma^2); %mode of distribution
%     headers={'P(GPa)','T(ºC)','t(hours)','meanDiameter','meanInterceptLength','modeFit','muFit','sigmaFit','EnstatiteFraction', 'NumberOfGrains'};
    data(1:numel(headers),ii)=[P_exp(ii); T_exp(ii); t_exp(ii); meanDiameter; meanInterceptLength; mod_d; mu; sigma; fEn; nGrains];
    disp([sampleNames{ii},': ',num2str(toc)])
end

    dlmwrite(filename,data,'delimiter','\t','precision',4,'-append')
    save('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/grainGrowthData.mat','data','headers','sampleNames')
%% Export .mat file to excel table
gg_data=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/grainGrowthData.mat');

sampleNames = {'FSG4','FSG5','Z1993','Z1962','Z1965','Z1968','B1272','B1273', 'B1274'...
    ,'B1275','B1276', 'Z2032', 'Z2033', 'Z2034', 'Z2035','Z2047', 'Z2049',...
    'Z2051','Z2060', 'A1178', 'A1179', 'Z2062','A1182', 'HV806', 'HV806-HP', 'HV806-HT'};

headers={'P_GPa','T_C','t_hours','d','dLIM','MoFit','uFit','sFit','fPx', 'n'};
gg_data=round(gg_data.data',2);

T = array2table(gg_data,'VariableNames',headers,'RowNames',sampleNames);
writetable(T,'/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/grainGrowthData.xls','WriteRowNames',true)

    %% Histograms
% histogram in log scale fitted to a log normal distribution
close all
clearvars
% colormaps
all_t=[8, 12, 24, 48, 72];
all_T =[1050, 1200, 1400, 1520];
all_P=[1, 5, 7, 10, 12];
all_Pxf=[4, 10];
cmap_t=cmap(length(all_t));
cmap_T=cmap(length(all_T));
cmap_P=cmap(length(all_P));
cmap_Pxf=cmap(length(all_Pxf));
alpha_l=0.6; % transparency for line plot

% paths
basePath_6 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/6pctEn';
basePath_13 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/13pctEn';

% Binning
min_x=-1; max_x=2;
nBins=100;
binEdge=10.^(min_x:(max_x-min_x)/nBins:max_x); % create bin edges with logarithmic scale  

% set figure and plots 
mS=100; %marker Size
fS=11; %Font size
lfS=9; %legend fontsize
phaseName='olivine';

%Uncertainties
P_err_Pc=0.02;% GPa

 %***************************----------------***************************
 %% Histograms starting material
 %Plot starting material
 f=figure('Units','normalized','pos',[0.25 0.25 0.2 0.25],'NumberTitle', 'off', 'Name', 'starting Material Histograms');
sp=1;
sampleNames={'FSG4','FSG5','Z1993', 'HV806'};
exp_fname={
                    fullfile(basePath_6,'startingMaterial')       % Starting material FSG4 (FSG4-GMF2)
                    fullfile(basePath_13,'startingMaterial')	  % Starting material FSG5
                    fullfile(basePath_6,'10gpa/0C-0h')           % Z1993 (MA3B) - No heating
%                     fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-StartingMaterial')    % HV806-Starting Material
                    };
color_sm=[lines(numel(exp_fname)),repmat(alpha_l,numel(exp_fname),1)];

    for ii=1: length(exp_fname)
        grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
        accumulatedDiameter=[]; %reset list for each sample
    
        for jj=1:length(grainsFileNameList)
            grains=load(grainsFileNameList{jj});
            grains=grains.grains;
            accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
        end

        pd= fitdist(accumulatedDiameter,'lognormal'); 
        pdf_counts = pdf(pd,binEdge);
        mu=pd.mu;sigma=pd.sigma;
        m_d=exp(mu+0.5*sigma^2); %mean
        md_d=exp(mu); % median
        mod_d=exp(mu-sigma^2); % mode
        
        hold on
        p=plot(binEdge,pdf_counts,'LineWidth',3,'Color',color_sm(ii,:),'DisplayName',sampleNames{ii});

        hold on
        scatter(m_d, pdf(pd,m_d),mS,'o','filled', 'MarkerEdgeColor', 'w','MarkerFaceColor',p.Color,'HandleVisibility','off');
        hold on
        scatter(md_d, pdf(pd,md_d),mS,'d', 'MarkerEdgeColor', 'w','MarkerFaceColor', p.Color,'HandleVisibility','off');
        hold on
        scatter(mod_d, pdf(pd,mod_d),mS,'s', 'MarkerEdgeColor','w','MarkerFaceColor',p.Color,'HandleVisibility','off');
        hold off
    end
makeScatterLegend
set(gca,'xscale','log'); % scale the x-axis
legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([0 30]); ylim([0 0.5]) 
axis square
hold off
box(gca,'on');

print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/Histogram_SM','-dpdf', '-painters')

%%  Grain growth data
% Time series
f=figure('Units','normalized','pos',[0.001 0.001 0.658 0.836],'NumberTitle', 'off', 'Name', 't, T & P histograms');
ids_scatterPlot=[3,3,3,6,6,6,9,9,9];
sp=1;

% Time series at 1GPa

exp_fname={   
                    fullfile(basePath_13,'1gpa/1400C-8h')      % A1182
                    fullfile(basePath_13,'1gpa/1400C-12h')    % Z2047
                    fullfile(basePath_13,'1gpa/1400C-24h')    % Z2060
                    fullfile(basePath_13,'1gpa/1400C-72h')    % Z2051           
                    };

t_series=[8; 12; 24; 72]; % Time in hours
P_series=[1; 1; 1; 1]; %Pressure in GPa

for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
              
    pd= fitdist(accumulatedDiameter,'lognormal'); 
    pdf_counts = pdf(pd,binEdge);

    subplot(3,3,sp);
    hold on
    p=plot(binEdge,pdf_counts,'LineWidth',3,'Color',[cmap_t(t_series(ii)==all_t,:),alpha_l],'DisplayName',[num2str(t_series(ii)),' h']);

    mu=pd.mu;sigma=pd.sigma;
    m_d=exp(mu+0.5*sigma^2);
    md_d=exp(mu);
    mod_d=exp(mu-sigma^2);

    hold on
    scatter(m_d, pdf(pd,m_d),mS,'o','filled', 'MarkerEdgeColor', 'w','MarkerFaceColor',cmap_t(t_series(ii)==all_t,:),'HandleVisibility','off','DisplayName',' Mean');
    hold on
    scatter(md_d, pdf(pd,md_d),mS,'d', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_t(t_series(ii)==all_t,:),'HandleVisibility','off','DisplayName',' Median');
    hold on
    scatter(mod_d, pdf(pd,mod_d),mS,'s', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_t(t_series(ii)==all_t,:),'HandleVisibility','off','DisplayName',' Mode');
    hold off
    
    %Scatter plot
    subplot(3,3,ids_scatterPlot(sp));
    hold on
    s=scatter(t_series(ii), mod_d,mS,'s', 'MarkerEdgeColor','k','DisplayName','1 GPa');
    hold off
   if ii>1
        s.HandleVisibility='off';
    end
    
end

subplot(3,3,sp);
makeScatterLegend
set(gca,'xscale','log'); % scale the x-axis
legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.4]) 
axis square
hold off

subplot(3,3,ids_scatterPlot(sp));
xlabel('Time (h)', 'FontSize',fS)
ylabel(' Grain size (µm)', 'FontSize',fS ) 
legend({}, 'Location','best','FontSize',lfS)
set(gca,'FontSize',lfS)
ylim([2 10]);
axis square
hold off

% Time series at 7GPa
exp_fname={
%                     fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_13,'7gpa/1400C-12h')    % Z2047
                    fullfile(basePath_13,'7gpa/1400C-24h')    % Z2060
                    fullfile(basePath_13,'7gpa/1400C-48h')    % Z2062
                    fullfile(basePath_13,'7gpa/1400C-72h')    % Z2051
                    };

sp=sp+1;
t_series=[12; 24; 48; 72]; % Time in hours
P_series=[7; 7; 7; 7]; %Pressure in GPa

for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
              
    pd= fitdist(accumulatedDiameter,'lognormal'); 
    pdf_counts = pdf(pd,binEdge);

    subplot(3,3,sp);
    hold on
    p=plot(binEdge,pdf_counts,'LineWidth',3,'Color',[cmap_t(t_series(ii)==all_t,:),alpha_l],'DisplayName',[num2str(t_series(ii)),' h']);

    mu=pd.mu;sigma=pd.sigma;
    m_d=exp(mu+0.5*sigma^2);
    md_d=exp(mu);
    mod_d=exp(mu-sigma^2);

    hold on
    scatter(m_d, pdf(pd,m_d),mS,'o','filled', 'MarkerEdgeColor', 'w','MarkerFaceColor',cmap_t(t_series(ii)==all_t,:),'HandleVisibility','off','DisplayName',' Mean');
    hold on
    scatter(md_d, pdf(pd,md_d),mS,'d', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_t(t_series(ii)==all_t,:),'HandleVisibility','off','DisplayName',' Median');
    hold on
    scatter(mod_d, pdf(pd,mod_d),mS,'s', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_t(t_series(ii)==all_t,:),'HandleVisibility','off','DisplayName',' Mode');
    hold off
    
    %Scatter plot
    subplot(3,3,ids_scatterPlot(sp));
    hold on
    s=scatter(t_series(ii), mod_d,mS,'s','filled','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerFaceAlpha',alpha_l,'MarkerEdgeAlpha',alpha_l, 'DisplayName','7 GPa');
    hold off
    if ii>1
        s.HandleVisibility='off';
    end
    
end
subplot(3,3,sp);
makeScatterLegend
set(gca,'xscale','log'); % scale the x-axis
legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.4]) 
axis square
hold off

subplot(3,3,ids_scatterPlot(sp));
xlabel('Time (h)', 'FontSize',fS)
ylabel(' Grain size (µm)', 'FontSize',fS ) 
legend({}, 'Location','best','FontSize',lfS)
set(gca,'FontSize',lfS)
ylim([2 10]);
axis square
hold off

%***************************----------------***************************

% Temperature series
sp=sp+1;% skip the scatter plot before next histogram
%Temperature series at 1GPa

exp_fname={
%                     fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_13,'1gpa/1050C-24h')    % 
                    fullfile(basePath_13,'1gpa/1200C-24h')    % 
                    fullfile(basePath_13,'1gpa/1400C-24h')    % 
                    };
T_series=[1050, 1200, 1400]; % Temperature in ºC
P_series=[1; 1; 1]; %Pressure in GPa
sp=sp+1;
%
for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
              
    pd= fitdist(accumulatedDiameter,'lognormal'); 
    pdf_counts = pdf(pd,binEdge);

    subplot(3,3,sp);
    hold on
    p=plot(binEdge, pdf_counts,'LineWidth',3,'Color',[cmap_T(T_series(ii)==all_T,:),alpha_l],'DisplayName',[num2str(T_series(ii)),' ºC']);
    
    mu=pd.mu;sigma=pd.sigma;
    m_d=exp(mu+0.5*sigma^2);
    md_d=exp(mu);
    mod_d=exp(mu-sigma^2);

    hold on
    scatter(m_d, pdf(pd,m_d),mS,'o','filled', 'MarkerEdgeColor', 'w','MarkerFaceColor',cmap_T(T_series(ii)==all_T,:),'HandleVisibility','off','DisplayName',' Mean');
    hold on
    scatter(md_d, pdf(pd,md_d),mS,'d', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_T(T_series(ii)==all_T,:),'HandleVisibility','off','DisplayName',' Median');
    hold on
    scatter(mod_d, pdf(pd,mod_d),mS,'s', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_T(T_series(ii)==all_T,:),'HandleVisibility','off','DisplayName',' Mode');
    hold off
    
    %Scatter plot
    subplot(3,3,ids_scatterPlot(sp));
    hold on
    s=scatter(T_series(ii), mod_d,mS,'s','MarkerEdgeColor','k','DisplayName',[num2str(P_series(ii)),' GPa']);
    hold off
    if ii>1
        s.HandleVisibility='off';
    end
    
end
subplot(3,3,sp);
makeScatterLegend
set(gca,'xscale','log'); % scale the x-axis
legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.4]) 
axis square

subplot(3,3,ids_scatterPlot(sp));
xlabel('Temperature (ºC)', 'FontSize',fS)
ylabel(' Grain size (µm)', 'FontSize',fS ) 
legend({}, 'Location','best','FontSize',lfS)
set(gca,'FontSize',lfS)
ylim([2 10]);
axis square
% Temperature series at 7GPa

exp_fname={
%                     fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_13,'7gpa/1050C-24h')    % 
                    fullfile(basePath_13,'7gpa/1200C-24h')    % 
                    fullfile(basePath_13,'7gpa/1400C-24h')    % 
                    fullfile(basePath_13,'7gpa/1520C-24h')    % 
                    };

sp=sp+1;
T_series=[1050, 1200, 1400, 1520]; % Temperature in ºC
P_series=[7; 7; 7; 7]; %Pressure in GPa

for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
              
    pd= fitdist(accumulatedDiameter,'lognormal'); 
    pdf_counts = pdf(pd,binEdge);

    subplot(3,3,sp);
    hold on
    p=plot(binEdge,pdf_counts,'LineWidth',3,'Color',[cmap_T(T_series(ii)==all_T,:),alpha_l],'DisplayName',[num2str(T_series(ii)),' ºC']);

    mu=pd.mu;sigma=pd.sigma;
    m_d=exp(mu+0.5*sigma^2);
    md_d=exp(mu);
    mod_d=exp(mu-sigma^2);

    hold on
    scatter(m_d, pdf(pd,m_d),mS,'o','filled', 'MarkerEdgeColor', 'w','MarkerFaceColor',cmap_T(T_series(ii)==all_T,:),'HandleVisibility','off','DisplayName',' Mean');
    hold on
    scatter(md_d, pdf(pd,md_d),mS,'d', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_T(T_series(ii)==all_T,:),'HandleVisibility','off','DisplayName',' Median');
    hold on
    scatter(mod_d, pdf(pd,mod_d),mS,'s', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_T(T_series(ii)==all_T,:),'HandleVisibility','off','DisplayName',' Mode');
    hold off
    
    %Scatter plot
    subplot(3,3,ids_scatterPlot(sp));
    hold on
    s=scatter(T_series(ii), mod_d,mS,'s','filled', 'MarkerFaceColor','k','MarkerEdgeColor','w','MarkerFaceAlpha',alpha_l,'MarkerEdgeAlpha',alpha_l,'DisplayName',[num2str(P_series(ii)),' GPa']);
    hold off
    if ii>1
        s.HandleVisibility='off';
    end
    
end
subplot(3,3,sp);
makeScatterLegend
set(gca,'xscale','log'); % scale the x-axis
legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.4]) 
axis square
hold off

subplot(3,3,ids_scatterPlot(sp));
xlabel('Temperature (ºC)', 'FontSize',fS)
ylabel(' Grain size (µm)', 'FontSize',fS ) 
legend({}, 'Location','best','FontSize',lfS)
set(gca,'FontSize',lfS)
ylim([2 10]);
axis square
hold off

% Pressure series
sp=sp+1;% skip the scatter plot before next histogram

% at 1400ºC, 24h, 4% En
exp_fname={
%                     fullfile(basePath_6,'startingMaterial')    % Starting material FSG4
                    fullfile(basePath_6,'5gpa/1400C-24h')   % 
                    fullfile(basePath_6,'7gpa/1400C-24h')   % 
                    fullfile(basePath_6,'12gpa/1400C-24h')    % 
                    };
P_series=[5; 7; 12];  % Pressure in GPa
Pxf_series=[4; 4; 4];

sp=sp+1;

for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
              
    pd= fitdist(accumulatedDiameter,'lognormal'); 
    pdf_counts = pdf(pd,binEdge);

    subplot(3,3,sp);
    hold on
    p=plot(binEdge,pdf_counts,'LineWidth',3,'Color',[cmap_P(P_series(ii)==all_P,:),alpha_l],'DisplayName',[num2str(P_series(ii)),' GPa']);

    mu=pd.mu;sigma=pd.sigma;
    m_d=exp(mu+0.5*sigma^2);
    md_d=exp(mu);
    mod_d=exp(mu-sigma^2);

    hold on
    scatter(m_d, pdf(pd,m_d),mS,'o','filled', 'MarkerEdgeColor', 'w','MarkerFaceColor',cmap_P(P_series(ii)==all_P,:),'HandleVisibility','off','DisplayName',' Mean');
    hold on
    scatter(md_d, pdf(pd,md_d),mS,'d', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_P(P_series(ii)==all_P,:),'HandleVisibility','off','DisplayName',' Median');
    hold on
    scatter(mod_d, pdf(pd,mod_d),mS,'s', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_P(P_series(ii)==all_P,:),'HandleVisibility','off','DisplayName',' Mode');
    hold off
    
    %Scatter plot
    subplot(3,3,ids_scatterPlot(sp));
    hold on
    s=scatter(P_series(ii), mod_d,mS,'s','MarkerEdgeColor','k','DisplayName',[num2str(Pxf_series(ii)),' vol.% Px']);
    hold off
    if  ii>1
        s.HandleVisibility='off';
    end
    
end
subplot(3,3,sp);
makeScatterLegend
set(gca,'xscale','log'); % scale the x-axis
legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.4]) 
axis square
hold off

subplot(3,3,ids_scatterPlot(sp));
xlabel('Pressure (GPa)', 'FontSize',fS)
ylabel(' Grain size (µm)', 'FontSize',fS ) 
legend({}, 'Location','best','FontSize',lfS)
set(gca,'FontSize',lfS)
ylim([2 10]);
axis square
hold off

% at 1400ºC, 24h, 10% En
exp_fname={
%                     fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_13,'1gpa/1400C-24h')   % B1272
                    fullfile(basePath_13,'7gpa/1400C-24h')   % Z2049
                    fullfile(basePath_13,'10gpa/1400C-24h')    % Z2032
                    };
P_series=[1; 7; 10];  % Pressure in GPa
Pxf_series=[10; 10; 10];
sp=sp+1;
subplot(3,3,sp);

for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
              
    pd= fitdist(accumulatedDiameter,'lognormal'); 
    pdf_counts = pdf(pd,binEdge);

    subplot(3,3,sp);
    hold on
    p=plot(binEdge,pdf_counts,'LineWidth',3,'Color',[cmap_P(P_series(ii)==all_P,:),alpha_l],'DisplayName',[num2str(P_series(ii)),' GPa']);

    mu=pd.mu;sigma=pd.sigma;
    m_d=exp(mu+0.5*sigma^2);
    md_d=exp(mu);
    mod_d=exp(mu-sigma^2);

    hold on
    scatter(m_d, pdf(pd,m_d),mS,'o','filled', 'MarkerEdgeColor', 'w','MarkerFaceColor',cmap_P(P_series(ii)==all_P,:),'HandleVisibility','off','DisplayName',' Mean');
    hold on
    scatter(md_d, pdf(pd,md_d),mS,'d', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_P(P_series(ii)==all_P,:),'HandleVisibility','off','DisplayName',' Median');
    hold on
    scatter(mod_d, pdf(pd,mod_d),mS,'s', 'MarkerEdgeColor', 'w','MarkerFaceColor', cmap_P(P_series(ii)==all_P,:),'HandleVisibility','off','DisplayName',' Mode');
    hold off
    
    %Scatter plot
    subplot(3,3,ids_scatterPlot(sp));
    hold on
    s=scatter(P_series(ii), mod_d,mS,'s','filled','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerFaceAlpha',alpha_l,'MarkerEdgeAlpha',alpha_l, 'DisplayName',[num2str(Pxf_series(ii)),' vol.% Px']);
    hold off
    if ii>1
        s.HandleVisibility='off';
    end
    
end
subplot(3,3,sp);
makeScatterLegend
set(gca,'xscale','log'); % scale the x-axis
legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.4]) 
axis square
hold off

subplot(3,3,ids_scatterPlot(sp));
xlabel('Pressure (GPa)', 'FontSize',fS)
ylabel(' Grain size (µm)', 'FontSize',fS ) 
legend({}, 'Location','best','FontSize',lfS)
set(gca,'FontSize',lfS)
ylim([2 10]);
axis square
hold off

lettersIds=fliplr({'a)','c)','b)','d)','f)','e)','g)','i)','h)'});
axs = findobj(f, 'type', 'axes', '-or', 'type', 'matlab.graphics.axis.Axes');

for axxN=1:length(lettersIds)
    box(axs(axxN),'on');
    t=title(axs(axxN),lettersIds(axxN),'FontSize',13);
    set(t, 'horizontalAlignment', 'left')
    set(t, 'units', 'normalized')
    t.Position(1)=0;
end

print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/Histograms','-dpdf', '-painters')

%% Histogram showing only data
clearvars
close all
phaseName='olivine';
mS=100; %marker Size
fS=11; %Font size
lfS=9;
basePath_6 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/6pctEn';
basePath_13 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/13pctEn';
gg_data=load('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/grainGrowthData.mat');
headers=gg_data.headers;
gg_data=gg_data.data;

P_exp=gg_data(1,:);
T_exp=gg_data(2,:);
T_exp(15)=1520;
t_exp=gg_data(3,:);
meanDiameter=gg_data(4,:);
meanIntercept=gg_data(5,:);
meanfEn=gg_data(6,:);
nGrains=gg_data(7,:);
f=figure('Units','normalized','pos',[0.001 0.001 0.693 0.98],'NumberTitle', 'off', 'Name', 'wholePage');

% Histograms
% histogram in log scale fitted to a log normal distribution

min_x=-1; max_x=2;
nBins=100;
binEdge=10.^(min_x:(max_x-min_x)/nBins:max_x); % create bin edges with logarithmic scale  
sp=1;
subplot(3,2,sp);

%***************************----------------***************************

% Time series at 1GPa
exp_fname={
                    fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_13,'1gpa/1400C-12h')    % Z2047
                    fullfile(basePath_13,'1gpa/1400C-24h')    % Z2060
                    fullfile(basePath_13,'1gpa/1400C-72h')    % Z2051
                    };
cmapBlue2Red=cmap(length(exp_fname));

t_series=[0; 12; 24; 72]; % Time in hours


for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
        
    hold on
    h=histogram(accumulatedDiameter,binEdge,'Normalization','pdf','FaceColor',cmapBlue2Red(ii,:),'EdgeColor',cmapBlue2Red(ii,:),'FaceAlpha',0.10,'EdgeAlpha',0.2,'DisplayName',[num2str(t_series(ii)), ' h']); % create the plot

    if ii==1
    %             h.DisplayName='Starting Material';
            h.DisplayName= 'Starting Material';
    end
        
end

set(gca,'xscale','log'); % scale the x-axis

legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.4]) 
axis square

% Time series at 7GPa

exp_fname={
                    fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_13,'7gpa/1400C-12h')    % Z2047
                    fullfile(basePath_13,'7gpa/1400C-24h')    % Z2060
                    fullfile(basePath_13,'7gpa/1400C-72h')    % Z2051
                   };

sp=sp+1;
subplot(3,2,sp);

for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
        
    hold on
    h=histogram(accumulatedDiameter,binEdge,'Normalization','pdf','FaceColor',cmapBlue2Red(ii,:),'EdgeColor',cmapBlue2Red(ii,:),'FaceAlpha',0.10,'EdgeAlpha',0.2,'DisplayName',[num2str(t_series(ii)), ' h']); % create the plot

    if ii==1
    %             h.DisplayName='Starting Material';
            h.DisplayName= 'Starting Material';
    end

end
    set(gca,'xscale','log'); % scale the x-axis

legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize',fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.4]) 
axis square

%***************************----------------***************************

% Temperature series
%Temperature series at 1GPa
T_series=[0, 1050, 1200, 1400, 1520]; % Temperature in ºC
cmapBlue2Red=cmap(length(T_series));

exp_fname={
                    fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_13,'1gpa/1050C-24h')    % 
                    fullfile(basePath_13,'1gpa/1200C-24h')    % 
                    fullfile(basePath_13,'1gpa/1400C-24h')    % 
                    };
sp=sp+1;
subplot(3,2,sp);

for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
        
    hold on
    h=histogram(accumulatedDiameter,binEdge,'Normalization','pdf','FaceColor',cmapBlue2Red(ii,:),'EdgeColor',cmapBlue2Red(ii,:),'FaceAlpha',0.10,'EdgeAlpha',0.2,'DisplayName',[num2str(T_series(ii)), ' ºC']); % create the plot

    if ii==1
            h.DisplayName= 'Starting Material';
    end

end

    set(gca,'xscale','log'); % scale the x-axis

legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize',fS)
set(gca,'FontSize',lfS)

xlim([1 30]); ylim([0 0.4])
axis square


% Temperature series at 7GPa

exp_fname={
                fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                fullfile(basePath_13,'7gpa/1050C-24h')    % Z2047
                fullfile(basePath_13,'7gpa/1200C-24h')    % Z2060
                fullfile(basePath_13,'7gpa/1400C-24h')    % Z2051
                fullfile(basePath_13,'7gpa/1520C-24h')    % Z2051
                };

sp=sp+1;
subplot(3,2,sp);
    
for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
  hold on
    h=histogram(accumulatedDiameter,binEdge,'Normalization','pdf','FaceColor',cmapBlue2Red(ii,:),'EdgeColor',cmapBlue2Red(ii,:),'FaceAlpha',0.10,'EdgeAlpha',0.2,'DisplayName',[num2str(T_series(ii)), ' ºC']); % create the plot

    if ii==1
            h.DisplayName= 'Starting Material';
    end
end

set(gca,'xscale','log'); % scale the x-axis

legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize',fS)
set(gca,'FontSize',lfS)

xlim([1 30]); ylim([0 0.4])  
axis square

%***************************----------------***************************

% Pressure series

% at 1400ºC, 24h
exp_fname={
                    fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_13,'1gpa/1400C-24h')   % B1272
                    fullfile(basePath_13,'7gpa/1400C-24h')   % Z2049
                    fullfile(basePath_13,'10gpa/1400C-24h')    % Z2032
                    };
P_series=[0; 1; 7; 10];  % Pressure in GPa
cmapBlue2Red=cmap(length(P_series));

sp=sp+1;
subplot(3,2,sp);
    
for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
        
    hold on
    h=histogram(accumulatedDiameter,binEdge,'Normalization','pdf','FaceColor',cmapBlue2Red(ii,:),'EdgeColor',cmapBlue2Red(ii,:),'FaceAlpha',0.10,'EdgeAlpha',0.2,'DisplayName',[num2str(P_series(ii)), ' GPa']); % create the plot

    if ii==1
            h.DisplayName= 'Starting Material';
    end
    
end
hold on

    set(gca,'xscale','log'); % scale the x-axis

legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize',fS)
set(gca,'FontSize',lfS)

xlim([1 30]); ylim([0 0.4])  
axis square

%***************************----------------***************************

% En content series
% at 1400ºC, 24h
exp_fname={
%                     fullfile(basePath_13,'startingMaterial')    % Starting material FSG5
                    fullfile(basePath_6,'7gpa/1400C-24h')    % MA2
                    fullfile(basePath_13,'7gpa/1400C-24h')   % Z2049
                    };
EnContent=[4,10];  % Pressure in GPa
cmapBlue2Red=cmap(length(EnContent));

sp=sp+1;
subplot(3,2,sp);
    
for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    accumulatedDiameter=[]; %reset list for each sample
    for jj=1:length(grainsFileNameList)
        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
        accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
    end
        
    hold on
    h=histogram(accumulatedDiameter, binEdge,'Normalization','pdf','FaceColor',cmapBlue2Red(ii,:),'EdgeColor',cmapBlue2Red(ii,:),'FaceAlpha',0.10,'EdgeAlpha',0.2,'DisplayName',[num2str(EnContent(ii)), ' vol.% Px']); % create the plot
    
end
hold on

    set(gca,'xscale','log'); % scale the x-axis

legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize',fS)
set(gca,'FontSize',lfS)

xlim([1 30]); ylim([0 0.4])  
axis square


lettersIds=fliplr({'a)','b)','c)','d)','e)','f)'});
axs = findobj(f, 'type', 'axes', '-or', 'type', 'matlab.graphics.axis.Axes');

for axxN=1:length(lettersIds)
    box(axs(axxN),'on');
    title(axs(axxN),lettersIds(axxN),'FontSize',13);
end

print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/Histograms_onlyData','-dpdf', '-painters')

%% Phase and IPF maps 
ol_color=[97, 148, 63]./255;
en_color=[98, 63, 35]./255;
nIdx_color = [1, 1, 1];
setMTEXpref('figsize','small')
sizeSquare=35; %um
% paths
basePath_6 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/6pctEn';
basePath_13 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/13pctEn';

exp_fname={

    fullfile(basePath_6,'startingMaterial')       % Starting material FSG4 (FSG4-GMF2)
    fullfile(basePath_13,'startingMaterial')	  % Starting material FSG5
    fullfile(basePath_6,'10gpa/0C-0h')           % Z1993 (MA3B) - No heating

    fullfile(basePath_6,'5gpa/1400C-24h')      % Z1962
    fullfile(basePath_6,'7gpa/1400C-24h')      % Z1965
    fullfile(basePath_6,'12gpa/1400C-24h')    % Z1968
    
    fullfile(basePath_13,'1gpa/1400C-24h')    % B1272
    fullfile(basePath_13,'1gpa/1200C-72h')    % B1273
    fullfile(basePath_13,'1gpa/1200C-24h')    % B1274
    
    fullfile(basePath_13,'1gpa/1050C-24h')    % B1275
    fullfile(basePath_13,'1gpa/1200C-12h')    % B1276
    fullfile(basePath_13,'10gpa/1400C-24h')  % Z2032
    
    fullfile(basePath_13,'7gpa/1200C-24h')    % Z2033
    fullfile(basePath_13,'7gpa/1200C-12h')    % Z2034
    fullfile(basePath_13,'7gpa/1050C-24h')    % Z2035
    
    fullfile(basePath_13,'7gpa/1400C-12h')    % Z2047
    fullfile(basePath_13,'7gpa/1520C-24h')    % Z2049
    fullfile(basePath_13,'7gpa/1400C-72h')    % Z2051
    
    fullfile(basePath_13,'7gpa/1400C-24h')    % Z2060
    fullfile(basePath_13,'1gpa/1400C-12h')    % A1178
    fullfile(basePath_13,'1gpa/1400C-72h')    % A1179
    
    fullfile(basePath_13,'7gpa/1400C-48h')    % Z2062
    fullfile(basePath_13,'1gpa/1400C-8h')      % A1182
    
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-StartingMaterial')    % HV806-Starting Material
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-HighPressure')      % HV806- High Pressure
    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-HighTemperature')    % HV806- High Temperature
};

sampleNames = {'FSG4-SM','FSG5-SM','Z1993','Z1962','Z1965','Z1968','B1272','B1273', 'B1274'...
    ,'B1275','B1276', 'Z2032', 'Z2033', 'Z2034', 'Z2035','Z2047', 'Z2049',...
    'Z2051','Z2060', 'A1178', 'A1179', 'Z2062','A1182', 'HV806-SM', 'HV806-HP', 'HV806-HT'};


for ii=1:length(exp_fname)
    grainsFileNameList= listEBSD(exp_fname{ii},'.mat');
    ebsdFileNameList= listEBSD(exp_fname{ii},'.ang');

    for jj=1:length(grainsFileNameList)
        [fPt,sampleName,~] = fileparts(grainsFileNameList{jj});

        grains=load(grainsFileNameList{jj});
        grains=grains.grains;
         
        grains('ol').color=ol_color;
        try 
            grains('en').color=en_color;
            grains('not').color=nIdx_color;
        catch
        end
        
        f=newMtexFigure;
        plot(grains,'micronbar','off')
        try
            g_ctd=grains.centroid;
            mean_x=mean(g_ctd(:,1));
            mean_y=mean(g_ctd(:,2));

            axis([mean_x-sizeSquare/2, mean_x+sizeSquare/2, mean_y-sizeSquare/2, mean_y+sizeSquare/2])
            setMTEXpref('xAxisDirection','east');setMTEXpref('zAxisDirection','intoPlane'); %EDAX: A1 up(-y), A2 left(-x)
            legend('off')
        catch
        end
        saveFigure(fullfile(fPt,['smallPhase_',sampleName,'_35x35um.png']))
    end
    
    for jj=1:length(ebsdFileNameList)
        [fPt,sampleName,~] = fileparts(ebsdFileNameList{jj});

        ebsd = EBSD.load(ebsdFileNameList{jj},'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');
          
        f=newMtexFigure;
        plot(ebsd,ebsd.prop.iq,'micronbar','off')
        mtexColorMap gray
        mean_x=mean(ebsd.prop.x(:));
        mean_y=mean(ebsd.prop.y(:));

        axis([mean_x-sizeSquare/2, mean_x+sizeSquare/2, mean_y-sizeSquare/2, mean_y+sizeSquare/2])
        setMTEXpref('xAxisDirection','east');setMTEXpref('zAxisDirection','intoPlane'); %EDAX: A1 up(-y), A2 left(-x)
        legend('off')
           
        saveFigure(fullfile(fPt,['small_iq_',sampleName,'_35x35um.png']))
    end
    
end
%% Show map for cleaning steps
close all
fname_raw='/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/ebsdCleaningExample/map20200201113816871_raw.ang';% raw data
fname_rescan='/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/ebsdCleaningExample/map20200201113816871_Rescan.ang'; %after re-scan
fname_clean_ps='/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/ebsdCleaningExample/map20200201113816871_Rescan_clean_ps.ang';% after rescan + dilation+ ps correction

fnames={fname_raw,fname_rescan,fname_clean_ps};

cs=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');
oM = ipfHSVKey(cs);
oM.inversePoleFigureDirection = zvector;% 001 IPF colorcoding 
phaseName='olivine';
f=newMtexFigure;
f.nrows=1;f.ncols=3;
setMTEXpref('figSize','small');
f.outerPlotSpacing=20;
lettersIds={'a)','b)','c)','d)','e)','f)','g)','h)','i)'};

for ii=1:length(fnames)
    fname=fnames{ii};
    ebsd =EBSD.load(fname,'interface','ang','convertEuler2SpatialReferenceFrame','setting 2');
    phaseNumber=find(strcmp(ebsd.mineralList,'olivine'))-1;

     if ii==length(fnames)
         [ebsd_clean,grains_clean]=ebsdCleaner(ebsd);
        ebsd_phase=ebsd_clean(ebsd_clean.phase==phaseNumber);
        plot(ebsd,ebsd.prop.iq,'micronbar','off')
        mtexColorMap gray
        colors = oM.orientation2color(ebsd_phase.orientations);
        hold on; plot(ebsd_phase, colors, 'FaceAlpha',0.7);% plot orientations
        hold on; plot(grains_clean('indexed').boundary,'k', 'linewidth', 0.9, 'edgealpha', 0.6);hold off;
     else
        ebsd_phase=ebsd(ebsd.phase==phaseNumber);
        plot(ebsd,ebsd.prop.iq,'micronbar','off')
        mtexColorMap gray
        colors = oM.orientation2color(ebsd_phase.orientations);
        hold on; plot(ebsd_phase, colors, 'FaceAlpha',0.7);% plot orientations
     end
     mtexTitle(lettersIds(ii),'alignLeftOutside', 'FontSize',14)
     hold off
    if ii~=length(fnames)
        f.nextAxis;
    end
end

saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/ebsdCleaningExample/cleaning.png')
%% Histograms HV806
close all
colors=lines(3);
 f=figure('Units','normalized','pos',[0.25 0.25 0.219 0.279],'NumberTitle', 'off', 'Name', 'HV806 Histograms');
sp=1;
sampleNames={'HV806','HV806-HT', 'HV806-HP'};
exp_fname={ 
                    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-StartingMaterial')    % HV806-Starting Material
                    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-HighTemperature')    % HV806-Starting Material
                    fullfile('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/HV806/HV806-HighPressure')    % HV806-Starting Material
                    };
    for ii=1: length(exp_fname)
        grainsFileNameList= listEBSD(exp_fname{ii},'.mat')
        accumulatedDiameter=[]; %reset list for each sample
    
        for jj=1:length(grainsFileNameList)
            grains=load(grainsFileNameList{jj});
            grains=grains.grains;
            accumulatedDiameter=[accumulatedDiameter; grains(phaseName).diameter];
        end
        length(accumulatedDiameter)
        pd= fitdist(accumulatedDiameter,'lognormal');
        pdf_counts = pdf(pd,binEdge);
        mu=pd.mu;sigma=pd.sigma;
        m_d=exp(mu+0.5*sigma^2); %mean
        md_d=exp(mu); % median
        mod_d=exp(mu-sigma^2) % mode
        hold on
        h=histogram(accumulatedDiameter, binEdge,'Normalization','pdf','FaceColor',colors(ii,:),'FaceAlpha',0.30,'EdgeAlpha',0.3,'DisplayName',sampleNames{ii}); % create the plot

        hold on
        p=plot(binEdge,pdf_counts,'LineWidth',3,'Color',colors(ii,:),'HandleVisibility','off');

        hold on
        scatter(m_d, pdf(pd,m_d),mS,'o','filled', 'MarkerEdgeColor', 'w','MarkerFaceColor',p.Color,'HandleVisibility','off');
        hold on
        scatter(md_d, pdf(pd,md_d),mS,'d', 'MarkerEdgeColor', 'w','MarkerFaceColor', p.Color,'HandleVisibility','off');
        hold on
        scatter(mod_d, pdf(pd,mod_d),mS,'s', 'MarkerEdgeColor','w','MarkerFaceColor',p.Color,'HandleVisibility','off');
        hold off
    end
makeScatterLegend
set(gca,'xscale','log'); % scale the x-axis
legend({},'FontSize',lfS,'Location','best')
setLogXticks
xlabel('Grain size  [μm]','FontSize',fS); ylabel('Probability Density','FontSize', fS)
set(gca,'FontSize',lfS)
xlim([1 30]); ylim([0 0.3]) 
axis square
hold off
print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/Histogram_HV806','-dpdf', '-painters')
%%
% hist grain size x neigbor ratio density plot
close all
f0=figure('Units','normalized','pos',[0.01 0.01 0.5 0.4],'NumberTitle', 'off', 'Name', 'Neighbor ratio density');

basePath_6 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/6pctEn';
basePath_13 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/13pctEn';

exp_fname={
%     fullfile(basePath_6,'5gpa/1400C-24h')      % Z1962
%     fullfile(basePath_13,'7gpa/1520C-24h')    % Z2049
%     fullfile(basePath_13,'1gpa/1400C-72h')    % A1179
       
%         fullfile(basePath_6,'startingMaterial')       % Starting material FSG4 (FSG4-GMF2)
        fullfile(basePath_13,'startingMaterial')	  % Starting material FSG5
    }
n_sample=1;% Sample every n grains

nBins=10;
minD=0;%0
maxD=9;%10
minNeighborRatio=0;
maxNeighborRatio=1;
N=[];edges=[];
x_binEdge=linspace(minD,maxD, nBins+1);
x_binCenter=x_binEdge(1:end-1)+diff(x_binEdge)/2;
y_binEdge=linspace(minNeighborRatio, maxNeighborRatio, nBins+1);
y_binCenter=y_binEdge(1:end-1)+diff(y_binEdge)/2;

for ii=1: length(exp_fname)
%     subplot(1,length(exp_fname),ii)
        grainsFileNameList= listEBSD(exp_fname{ii},'.mat')
        all_gsz=[]; all_neighborRatio=[];
        for jj=1:length(grainsFileNameList)
             
            grains=load(grainsFileNameList{jj});
            grains=grains.grains;
            disp('Min/ max grain size:')
            [min(grains.diameter), max(grains.diameter)]
            ol_phaseId=unique(grains('ol').phase);
            grains_id=grains.id;
            neighborRatio=nan(1,length(grains)); g_size=nan(1,length(grains));
            for g_id=1:n_sample:length(grains) % sample every n grains
                g_grainId=grains_id(g_id);
                g=grains('id', g_grainId);
                
                if g.phase==ol_phaseId
                    g_size(g_id)=g.diameter;
                    gb_gId=unique(g.boundary.grainId);
                    gb_gId=gb_gId(ismember(gb_gId,grains_id));
                    g_neigbr=grains('id',gb_gId(gb_gId~=0 & gb_gId~=g_grainId));
                    neighborRatio(g_id) = length(g_neigbr('En'))/length(g_neigbr);
                else
                    g_size(g_id)=nan;
                    neighborRatio(g_id)=nan;
                end
            end
            all_gsz=[all_gsz, g_size(~isnan(g_size))];
            all_neighborRatio=[all_neighborRatio, neighborRatio(~isnan(g_size))];  
        end
        % 2D density plot
        X = all_gsz';
        Y=all_neighborRatio';
        [mean(X), mean(Y)]
        for kk =1:nBins
            id_lines=(all_gsz>x_binEdge(kk) & all_gsz<x_binEdge(kk+1));
            h= histogram2(X(id_lines),Y(id_lines),'NumBins',[1 nBins],...
            'XBinEdges',[x_binEdge(kk) ,x_binEdge(kk+1)],'YBinEdges',y_binEdge,'Normalization','pdf');
            N(kk,:)=h.Values;
        end

        [xi,yi] = meshgrid(x_binEdge, y_binEdge);
        C=nan(size(N,1)+1,size(N,2)+1);
        C(1:end-1,1:end-1)=N';
        axs(ii)=pcolor(xi,yi,C);
        axis square tight
        set(gca,'colorscale','log')
        set(gca,'FontSize',14)

        xlabel('Grain Size [µm]','FontSize',16); ylabel('Neighbour Ratio','FontSize',16);
        axis square tight
        cb=colorbar;
        ylabel(cb,'Probability Density','FontSize',16);

end
% saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/neighbourRatio_FSG5SM.png');   
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/neighbourRatio_FSG5SM.pdf');     

%% Distance between pyroxene grains

close all
clearvars
f0=figure('Units','normalized','pos',[0.01 0.01 0.5 0.4],'NumberTitle', 'off', 'Name', 'Px distance');
basePath_6 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/6pctEn';
basePath_13 = '/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/ggEXP/13pctEn';

exp_fname={
%     fullfile(basePath_6,'5gpa/1400C-24h')      % Z1962
%     fullfile(basePath_13,'7gpa/1520C-24h')    % Z2049
%     fullfile(basePath_13,'1gpa/1400C-72h')    % A1179
       
%         fullfile(basePath_6,'startingMaterial')       % Starting material FSG4 (FSG4-GMF2)
        fullfile(basePath_13,'1gpa/1400C-24h')
        fullfile(basePath_13,'7gpa/1400C-24h')	  
        fullfile(basePath_13,'10gpa/1400C-24h')	  

    }
all_P=[1, 5, 7, 10, 12];
P_series=[1,7,10];
cmap_P=cmap(length(all_P));

for ii=1: length(exp_fname)
%     subplot(1,length(exp_fname),ii)
       grainsFileNameList= listEBSD(exp_fname{ii},'.mat')
       PxDist=[]; nPxGrains=[]; totalArea=[];
        for jj=1:length(grainsFileNameList)
                         
            grains=load(grainsFileNameList{jj});
            grains=grains.grains;
            Px_ctd=grains('En').centroid; %center of grains
%             figure;
%             scatter(Px_ctd(:,1),Px_ctd(:,2))
            D = tril(squareform(pdist(Px_ctd)));% distance of a to b is equal to b to a
            D(D==0)=nan; %remove distances from grain to itself
            min_D=nanmin(D)/mean(grains('ol').diameter);
            PxDist=[PxDist, min_D];
            nPxGrains(jj)=length(grains('En'));
            totalArea(jj)=sum(grains('En').area);
        end
        rationGrainsArea=sum(nPxGrains)/sum(totalArea)
        hold on
        h(jj)=histogram(PxDist,'Normalization','probability','BinWidth',0.1,'DisplayName',[num2str(P_series(ii)),' GPa'],'FaceColor',cmap_P(P_series(ii)==all_P,:));
        meanPxDist=nanmean(PxDist)
        binV=h(jj).Values;
        binC=h(jj).BinEdges(1:end-1)+diff(h(jj).BinEdges)/2;
        [~,id]=min(abs(binC-meanPxDist))
        hold on
        s(jj)=scatter3(meanPxDist,binV(id),2,90,'c','MarkerEdgeColor','w','MarkerFaceColor',cmap_P(P_series(ii)==all_P,:),'HandleVisibility','off');
%         nanmean(PxDist);
        alpha([h(jj),s(jj)],0.6);
end
hold on
scatter(nan, nan,100,'k', 'DisplayName',' Mean');
hold off
legend
xlabel('Normalized Distance','FontSize',16);
ylabel('Probability Density','FontSize',16);
set(gca,'FontSize',14)

axis square
xlim([0 3])
% saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/hist_minDistance.png')
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/hist_minDistance.pdf')

%% Plot FTIR
f0=figure('Units','normalized','pos',[0.001 0.17 0.3 0.3],'NumberTitle', 'off', 'Name', 'FTIR');
fname='/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/data/FTIR/FSG4_GMF2.TXT';
xRangeBgd=[2050, 4450];% cm^-1
%Get data
T = dlmread(fname);

x=T(:,1);
y=T(:,2);

indBgd= x>xRangeBgd(1) & x < xRangeBgd(2);

%%Remove background
y_out= msbackadj(x(indBgd),y(indBgd),'WindowSize',20);

plot(x(indBgd),y_out,'k','LineWidth',0.5,'DisplayName','Background subtracted data') %background subtracted
axis square 
xlim([2500, 4000])% cm^-1
set (gca, 'xdir', 'reverse' )
ax = gca;
xlabel('Wavenumber (cm^{-1})','FontSize',12);
ylabel('Absorbance','FontSize',12);
set(gca,'FontSize',10)
print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/FTIR','-dpdf', '-painters')
saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/FTIR.png')
%%
% Fit data to grain growth equation (3D plot)
clearvars
% close all
figure('Units','normalized','pos',[0.11 0.01 0.6 0.3],'NumberTitle', 'off', 'Name', 'Fit -grain growth data');
%Import data
gg_data=load('grainGrowthData.mat');
headers=gg_data.headers;
sampleNames=gg_data.sampleNames;
gg_data=gg_data.data;

P_exp=gg_data(1,:);% Pressure (GPa)
T_exp=gg_data(2,:); %Temperature ºC
t_exp=gg_data(3,:); %Time (h)
fEn=gg_data(9,:); %enstatite fraction
meanDiameter=gg_data(6,:); % grain size in um
d_0=meanDiameter(strcmp(sampleNames,'FSG5-SM')); %starting grain size

id=find(T_exp==1400 & fEn>0.07); %Get points for 1, 7, and 10 GPa (at 1400ºC, fPx>7%)

% Check if correct samples were chosen
% T = array2table([P_exp(id);T_exp(id);t_exp(id)],'VariableNames',sampleNames(id)','RowName',{'P' ,'T','t'}); 
% disp(T)

% parameters
%-------------------------

R   = 8.314; % Gas constant (J . mol^-1 . K^-1) 
d0  = 2.294*10^-6; % Initial grain size (m)
T   = 1673; %Temperature (K)
t   = t_exp(id)'*3600;% Time (s)
% t= 24 *3600;% Time (s)
P=P_exp(id)'*10^9; % Pressure in Pa
d=meanDiameter(id)'*10^-6; %grain size in m
% d = ( k0.*exp(-(Ea+Va*P)./(R*T)) .* t + d0.^p) .^(1/p);

% fit
%---------------------------

fo = fitoptions('Method','NonlinearLeastSquares', 'Lower', [1e-10, 1e-10, 200E3, 3], 'Upper', [1e-5, 1e-5, 700E3, 6],'Normalize','off','Algorithm','Trust-Region','Display','iter','TolFun',1e-13,'TolX',1e-10);

FitStr = ['(',num2str(d0),'^n + a*x * exp(-(c+ b*y)/',num2str(R*T),') ).^(1/n)'];

ft = fittype(FitStr,'options',fo,'coefficients',{'a','b','c','n'},'dependent',{'z'},'independent',{'x','y'});

[curve,~] = fit([t,P],d,ft,'Robust', 'LAR'); % curve fitting
sprintf('k_0= %.2d, V*= %.2d, E*=%.2f n= %.2f', curve.a, curve.b, curve.c, curve.n)
% Plot data and Fit

% % plot 3D data:
subplot(1,2,1)
p=plot(curve,[t,P],d);%
xlabel('Time [s] ');ylabel('Pressure [Pa]'); zlabel('Grain size [m]');
shading interp; alpha(p(1),0.9); colormap gray
axis square
%Error bar
y1=P+0.1*P; % + Error Pressure
y2=P-0.1*P; % - Error Pressure
z1=d+0.1*d;% + Error grain size
z2=d-0.1*d; % - Error grain size
colors=lines;
hold on 
plot3([t(:),t(:)]', [y1(:),y2(:)]', [d(:),d(:)]','Color',colors(1,:),'Linewidth',2) % Plot Error Bars Y
hold on 
plot3([t(:),t(:)]', [P(:),P(:)]', [z1(:),z2(:)]','Color',colors(1,:),'Linewidth',2) % Plot Error Bars Z
hold off
%
% 2D: at t=24h
subplot(1,2,2)
n=curve.n; %grain size exponent 
k0= curve.a; % Pre exponential term (m^n s^-1)
Va=curve.b; % Activation volume
Ea=curve.c; %Activation Energy

T = 1400+273.15;
d0  = 2.29*10^-6; % Initial grain size (m)
t  = 24*3600; % Time (s)
P = [1:12]*10^9; %Pressure in Pa
% Grain growth equation
k=k0*exp(-(Ea+P.*Va)/(R*T)); % Rate Constant  (P dependent)
GSt=((d0^n)+(k.*t)).^(1/n); % Final Grain Size (m)

id=find(T_exp==1400 & t_exp==24 & fEn>0.07); % Get points for 1, 7, and 10 GPa (at 1400ºC, 24h, fPx>7%)

hold on
%Errorbar
yneg = 0.1*meanDiameter(id);
ypos = 0.1*meanDiameter(id);
xneg = 0.1*P_exp(id);
xpos = 0.1*P_exp(id);
% scatter(P_exp(id), meanDiameter(id))
errorbar(P_exp(id), meanDiameter(id), yneg, ypos, xneg, xpos, 'o')

plot(P*10^-9, GSt*10^6,'k')

hold off
% xlim([0.95 12])

legend('Data','Fit')
hold off
xlabel('Pressure [GPa]')
ylabel('Grain size [µm]')
axis square
box on
% print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/grainGrowthFit2','-dpdf', '-painters')

%% Viscosity profiles

%% Plot gs at different depths, and expected viscosity profile - Single viscosity, constant Stress
clearvars
close all
line_smbl={'-','--',':','-.'};
line_colors=lines;
set(0, 'DefaultLineLineWidth', 2.5);
f=figure('Units','normalized','pos',[0.11 0.01 0.45 0.3],'NumberTitle', 'off', 'Name', 'Viscosity estimation with depth');

% close all

%  Parameters:
D=[200 :10: 410]; %Depth in km
P=0.0357*D-0.6096; %Pressure (P) in GPa as a function of depth (D) in km (PREM -Anderson)
% T = (-4.5153e-07)*D.^4 +(5.3209e-04)*D.^3 + (-0.2266*D.^2) + (40.6121*D) +(-1.0926e+03); %Temperature (T) in C as a function of depth (D) in km (PREM -Anderson)
T =0.5122*D+ 1277; %Temperature (T) in C as a function of depth (D) in km (Katsura et al 2004)

R=8.314; % Gas constant(J / mol* K)
GS0=1e-5 ; % Initial grain size (m)
t_s=[1 10]'*1e6*365*24*3600;% Time (1, 10 Ma in seconds)
P_Pa=P*10^9; %P in Pa
P_MPa=P*10^3; %P in MPa
T_K = T+273.15; %T  in K

% Rheological parameters

% disGBS Ohuchi (2015)- dry/wet olivine
% %-----------
%         A_gbs = 10^-4.89; % Pre exponential
%         n_gbs = 3; % Stress exponent
%         m_gbs = 1.1; % grain size exponent
%         Q_gbs = 4.23e5; % Activation energy (J/mol)
%         V_gbs = 17.6e-6; % Activation volume (m^3/mol)
%         r_gbs=1.25; % Exponent for water fugacity
%         wf=50; % Water fugacity  % fH20 in MPa or COH in ppm H/Si
%         % d=  grain size in m
%         % P = pressure in MPa
%         % sigma = stress in MPa
%         % T in K

% disGBS-  Hansen (JGR 2011)
% %-----------
A_gbs = 10^4.8; % Pre exponential
n_gbs = 2.9; % Stress exponent
m_gbs = 0.7; % grain size exponent
V_gbs = 18e-6; % as used in Hansen et al. (JGR 2011)
Q_gbs = 4.45e5 - 300e6*V_gbs; % Activation energy (J/mol), has to be corrected for using the activation volume
%         V_gbs = 3e-6; % Activation volume (m^3/mol) (Depths> 240km)

% d =  grain size in um
% P = pressure in Pa
% sigma = stress in MPa
% T in K

% Diffusion creep Hirth and Kohlstedt, 2003 (and Hansen et al., 2011)
%----------------
A_dif = 10^(7.6);  % in micron^md.MPa^-nd.s^-1
n_dif = 1.0;% Stress exponent
m_dif = 3.0;% grain size exponent
Q_dif = 3.75e5; % in J/mol
V_dif = 6.e-6;  % in m^3/mol
% d =  grain size in um
% P = pressure in Pa
% sigma = stress in MPa
% T in K

% Dislocation creep Hirth and Kohlstedt (2003) - dry olivine
A_dis = 1.1e5;  % Arrhenius constant in MPa^-nr.s^-1
n_dis = 3.5;    % stress exponent
Q_dis = 5.3e5;  % activation energy in J/mol
V_dis = 14e-6;  % activation volume in m^3.mol^-1 (e.g., Karato and Rubie 1997)
m_dis = 0;
% d =  grain size in um
% P = pressure in Pa
% sigma = stress in MPa
% T in K
% Def. at Constant stress
sigma=1;% Stress in MPa (constant stress)
sigma_Pa=sigma*10^6; % Stress in Pa

% Grain growth parameters (this work)
n=3.88; %grain size exponent
E=608*10^3;% Activation Energy(J/mol)
k0= 2.11e-7; % Pre exponential term (m^n s^-1)
Va=[0, 4.3e-6]; % Activation volume (m3/mol)

str= sprintf('d_0 = %.0f µm',GS0*1E6); % initial grain size 


% loop at different activation volumes
for jj=1:length(Va)
    V=Va(jj);
    % Grain growth equation
    k=k0*exp(-(E+P_Pa.*V)./(R.*T_K)); % Rate Constant  (P dependent)
    GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (m)
    GSt_um=GSt*1e6; % m to um
    
    % Calculate viscosity
    %Calculate strain rates for disGBS and diffusion creep
    
    % disGBS Ohuchi (2015)- dry/wet olivine
    %eps_gbs=A_gbs.*(sigma.^n_gbs).*(GSt.^-m_gbs).*(wf^r_gbs).* exp(-(Q_gbs+P_MPa.*V_gbs)./(R.*T_K));
    
    % disGBS Hansen (2011) - dry olivine
    eps_gbs=A_gbs.*(sigma.^n_gbs).*(GSt_um.^-m_gbs).* exp(-(Q_gbs+P_Pa.*V_gbs)./(R.*T_K));
    
    % Diffusion creep Hirth and Kohlstedt, 2003 (and Hansen et al., 2011)
    eps_dif=A_dif.*(sigma.^n_dif).*(GSt_um.^-m_dif).* exp(-(Q_dif+P_Pa.*V_dif)./(R.*T_K));
    
    % Dislocation creep Hirth and Kohlstedt, 2003
    eps_dis=A_dis.*(sigma.^n_dis).*(GSt_um.^-m_dis).* exp(-(Q_dis+P_Pa.*V_dis)./(R.*T_K));
    
    % Sum strain rates
    epsSum=eps_dif+eps_gbs+eps_dis;
    
    
    % compute deformational work rate
        %DefWork = sigma_Pa.*epsSum;
    
    % compute eq. grain size ?
    
    
    % Plots
    % a) Pressure and temperature profile with depth
    ax1=subplot(1,3,1);
    plot(T,D,'k')
    xlabel('Temperature [ºC]')
    ylabel ('Depth [km]')
    set (gca,'Ydir','reverse','box','on','LineWidth',2,'FontSize',14,'PlotBoxAspectRatio',[1 1.5 1])
    ylim([200 410])
    
    % b) Grain size profile with depth
    ax2=subplot(1,3,2);
    hold on
    plot(GSt_um(1,:),D, 'LineStyle', line_smbl{2}, 'Color',line_colors(jj,:), 'DisplayName', sprintf('1Ma, V^* =%.1e', V))
    hold on
    plot(GSt_um(2,:),D, 'LineStyle', line_smbl{3}, 'Color',line_colors(jj,:), 'DisplayName', sprintf('10Ma, V^* =%.1e', V))
    hold off
    xlabel('Grain size [µm]')
    %ylabel ('Depth [km]')
    axis tight
    set (gca,'Ydir','reverse','YTickLabel',[],'box','on','LineWidth',2,'FontSize',14,'PlotBoxAspectRatio',[1 1.5 1])
    xlim([0 4500])
    t1=text(0.98, 0.98, str,'Units','normalized','BackgroundColor',[1 1 1 0.5],'EdgeColor',[0 0 0 0.2],'Linewidth',1,'FontSize',10','VerticalAlignment','top','HorizontalAlignment','right','Margin',0.01);

    %legend;
    
    % c) Viscosity profile with depth
    ax3=subplot(1,3,3);
    yyaxis left
    hold on
    %1Ma
    plot(sigma_Pa./epsSum(1,:), D, 'LineStyle', line_smbl{2}, 'Color',line_colors(jj,:), 'DisplayName', sprintf('V^*= %.1e', V));
    hold on
    %10Ma
    plot(sigma_Pa./epsSum(2,:), D, 'LineStyle', line_smbl{3}, 'Color',line_colors(jj,:), 'DisplayName', sprintf('V^*= %.1e', V));
    xlabel('\eta [Pa s]')
    axis tight
    set(gca,'xscale','log','Ydir','reverse','ycolor','k','YTick',[],'Xlim',[1e19 1e23],'XTick',[1e19  1e21  1e23],'PlotBoxAspectRatio',[1 1.5 1]);
    t2=text(0.98, 0.98, str,'Units','normalized','BackgroundColor',[1 1 1 0.5],'EdgeColor',[0 0 0 0.2],'Linewidth',1,'FontSize',10','VerticalAlignment','top','HorizontalAlignment','right','Margin',0.01);

    yyaxis right
    ylim([P(1), P(end)])
    ylabel ('Pressure [GPa]')
%     legend('Orientation','vertical')
    set(gca,'Ydir','reverse','ycolor','k','LineWidth',2,'FontSize',14)
    
end
lettersIds=fliplr({'a)','b)','c)'});
axs = findobj(f, 'type', 'axes', '-or', 'type', 'matlab.graphics.axis.Axes');

for axxN=1:length(lettersIds)
    box(axs(axxN),'on');
    t=title(axs(axxN),lettersIds(axxN),'FontSize',14,'FontWeight','normal');
    set(t, 'horizontalAlignment', 'left')
    set(t, 'units', 'normalized')
    t.Position(1)=0;
end
print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/viscosityPofile','-dpdf', '-painters')

%% Effect of different starting grain sizes
clearvars
line_smbl={'-','--',':','-.'};
line_colors=lines;
set(0, 'DefaultLineLineWidth', 2.5);

f=figure('Units','normalized','pos',[1.25 1 0.285 1.05],'NumberTitle', 'off', 'Name', 'Viscosity estimation with depth - Effect of dif. starting grain size');

%  Parameters:
D=[200 :10: 410]; %Depth in km
P=0.0357*D-0.6096; %Pressure (P) in GPa as a function of depth (d) in km (PREM -Anderson)
% T = (-4.5153e-07)*D.^4 +(5.3209e-04)*D.^3 + (-0.2266*D.^2) + (40.6121*D) +(-1.0926e+03); %Temperature (T) in C as a function of depth (d) in km (PREM -Anderson)
T =0.5122*D+ 1277; %Temperature (T) in C as a function of depth (D) in km (Katsura et al 2004)

R=8.314; % Gas constant(J / mol* K)
%     GS0=10*1e-6 ; % Initial grain size (m) 
t_s=[1 10]'*1e6*365*24*3600;% Time (1 Ma in seconds)
P_Pa=P*10^9; %P in Pa
P_MPa=P*10^3; %P in MPa
T_K = T+273.15; %T  in K
% Constant stress
sigma=1;% Stress in MPa (constant stress)
sigma_Pa=sigma*10^6; % Stress in Pa

% Grain growth parameters (this work)
n=3.88; %grain size exponent 
E=608*10^3;% Activation Energy(J/mol)
k0= 2.11e-7; % Pre exponential term (m^n s^-1)
Va=[0, 4.3e-6]; % Activation volume (m3/mol)
d_initial = [10, 100, 1000]*1e-6; % Initial grain size (m) 

% Deformation parameters

% disGBS-  Hansen (JGR 2011)
% %-----------
A_gbs = 10^4.8; % Pre exponential
n_gbs = 2.9; % Stress exponent
m_gbs = 0.7; % grain size exponent
V_gbs = 18e-6; % as used in Hansen et al. (JGR 2011)
Q_gbs = 4.45e5 - 300e6*V_gbs; % Activation energy (J/mol), has to be corrected for using the activation volume
%         V_gbs = 3e-6; % Activation volume (m^3/mol) (Depths> 240km)
% d =  grain size in um
% P = pressure in Pa
% sigma = stress in MPa
% T in K

% Diffusion creep Hirth and Kohlstedt, 2003 (and Hansen et al., 2011)
%----------------
A_dif = 10^(7.6);  % in micron^md.MPa^-nd.s^-1 
n_dif = 1.0;% Stress exponent
m_dif = 3.0;% grain size exponent
Q_dif = 3.75e5; % in J/mol
V_dif = 6.e-6;  % in m^3/mol
% d =  grain size in um    
% P = pressure in Pa
% sigma = stress in MPa
% T in K

% Dislocation creep Hirth and Kohlstedt (2003) - dry olivine
A_dis = 1.1e5;  % Arrhenius constant in MPa^-nr.s^-1
n_dis = 3.5;    % stress exponent
Q_dis = 5.3e5;  % activation energy in J/mol
V_dis = 14e-6;  % activation volume in m^3.mol^-1 (e.g., Karato and Rubie 1997)
m_dis = 0;
% d =  grain size in um    
% P = pressure in Pa
% sigma = stress in MPa
% T in K

nlines=length(d_initial); % number of rows in plot (number of different initial grain sizes)
for kk=1:nlines
    GS0 = d_initial(kk);
    str= sprintf('d_0 = %.0f µm',GS0*1E6);

    % loop at different activation volumes
    for jj=1:length(Va)
        V=Va(jj);

        % Grain growth equation
        k=k0*exp(-(E+P_Pa.*V)./(R.*T_K)); % Rate Constant  (P dependent)
        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (m)
        GSt_um=GSt*1e6; % m to um


        % Calculate viscosity

        % disGBS Ohuchi (2015)- dry/wet olivine
        %eps_gbs=A_gbs*(sigma.^n_gbs)*(GSt.^-m_gbs)*(wf^r_gbs)* exp(-(Q_gbs+P_MPa.*V_gbs)/(R.*T_K));

        % disGBS Hansen (2011) - dry olivine
        eps_gbs=A_gbs.*(sigma.^n_gbs).*(GSt_um.^-m_gbs).* exp(-(Q_gbs+P_Pa.*V_gbs)./(R.*T_K));

        % Diffusion creep Hirth and Kohlstedt, 2003 (and Hansen et al., 2011)
        eps_dif=A_dif.*(sigma.^n_dif).*(GSt_um.^-m_dif).* exp(-(Q_dif+P_Pa.*V_dif)./(R.*T_K));

        % Dislocation creep Hirth and Kohlstedt, 2003
        eps_dis=A_dis.*(sigma.^n_dis).*(GSt_um.^-m_dis).* exp(-(Q_dis+P_Pa.*V_dis)./(R.*T_K));

        % Sum strain rates
        epsSum=eps_gbs+eps_dif+eps_dis;

        % Plots

        % i) Grain size profile with depth
        subplot(nlines,2,kk*2-1);
        hold on
        plot(GSt_um(1,:),D, 'LineStyle', line_smbl{2}, 'Color',line_colors(jj,:), 'DisplayName', sprintf('1Ma, V^* = %.1e', V))
        hold on
        plot(GSt_um(2,:),D, 'LineStyle', line_smbl{3}, 'Color',line_colors(jj,:), 'DisplayName', sprintf('10Ma, V^* = %.1e', V))
        hold off
        xlabel('Grain size [µm]')
        ylabel ('Depth [km]')

        axis tight
        set (gca,'Ydir','reverse','YTickLabel',[],'box','on','LineWidth',2,'FontSize',14,'PlotBoxAspectRatio',[1 1.5 1])
        ylim([200 410])
        xlim([0 4500])
        t1=text(0.98, 0.98, str,'Units','normalized','BackgroundColor',[1 1 1 0.5],'EdgeColor',[0 0 0 0.2],'Linewidth',1,'FontSize',10','VerticalAlignment','top','HorizontalAlignment','right','Margin',0.01);

        % ii) Viscosity profile with depth
        subplot(nlines,2,kk*2);
        yyaxis left
        hold on
        %1Ma
        plot(sigma_Pa./epsSum(1,:), D, 'LineStyle', line_smbl{2}, 'Color',line_colors(jj,:), 'DisplayName', sprintf('V^* = %.1e', V));
        hold on
        %10Ma
        plot(sigma_Pa./epsSum(2,:), D, 'LineStyle', line_smbl{3}, 'Color',line_colors(jj,:), 'DisplayName', sprintf('V^* = %.1e', V));
        xlabel('\eta [Pa s]')
        axis tight
        set(gca,'xscale','log','Ydir','reverse','ycolor','k','YTick',[],'Xlim',[1e19 1e23],'XTick',[1e19  1e21  1e23],'PlotBoxAspectRatio',[1 1.5 1],'box','on');

        yyaxis right
        ylim([P(1), P(end)])
        ylabel ('Pressure [GPa]')
        set(gca,'Ydir','reverse','ycolor','k','LineWidth',2,'FontSize',14,'PlotBoxAspectRatio',[1 1.5 1],'box','on')
        t2=text(0.98, 0.98, str,'Units','normalized','BackgroundColor',[1 1 1 0.5],'EdgeColor',[0 0 0 0.2],'Linewidth',1,'FontSize',10','VerticalAlignment','top','HorizontalAlignment','right','Margin',0.01);
    end
    
end
    lettersIds=fliplr({'a)','b)','c)','d)','e)','f)'});
    axs = findobj(f, 'type', 'axes', '-or', 'type', 'matlab.graphics.axis.Axes');

    for axxN=1:length(lettersIds)
        box(axs(axxN),'on');
        t=title(axs(axxN), lettersIds(axxN),'FontSize',12,'FontWeight','normal');
        set(t, 'horizontalAlignment', 'left')
        set(t, 'units', 'normalized')
        t.Position(1)=0;
    end
    print('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/viscosityPofile_d0','-dpdf', '-painters')
%% Compare our experimental data with other grain growth Laws:

% Grain growth law
% k=k0*exp(-(E)/(R*T_K)); % Rate Constant
% GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)

%% Plot grain growth laws and experimental data:
close all

f=figure('Units','normalized','Position', [-0.0396 1.1009 1.1010 1.1722]);
subplot (2,1,1)

% Constants
R=8.314;%*10^-3; % Gas constant(J . mol^-1 . K^-1) 
% T=[1050,1200,1400,1520]; % Temperature (ºC)
T=[1400];
P_all=[1,5,7,10,12]; %GPa
V=0.7e-6;% m^3/mol % Si diffusion-Bejina 1999)
    
t_min=1;% minimum time (hours)
t_max=3;% maximum time (hours)
t_resolution= 50;% number of points for each line
t=logspace(t_min,t_max, t_resolution); % %time points (hours)
t_s=t*3600; % Time(s)

GS0_um=2.5
GS0=GS0_um*1e-6 ; % Initial grain size (m) 

lgd={};

%symbols/colors
% symbols: pressure
%not filled 0>enstatite>6%, filled: 6%>enstatite>15%
%Colors: temperature
length_cmap = length(T);% 1050, 1200, 1400ºC
blue = [29, 173, 207]/255;
red = [249, 84, 28]/255;
alpha_value=0.6;
cmapBlue2Red = [linspace(blue(1),red(1),length_cmap)', linspace(blue(2),red(2),length_cmap)', linspace(blue(3),red(3),length_cmap)'];
cmapBlue2Red_alpha = [linspace(blue(1),red(1),length_cmap)', linspace(blue(2),red(2),length_cmap)', linspace(blue(3),red(3),length_cmap)',repmat(alpha_value,length_cmap,1)];

mk_smbl={'o','s','d','^','p'}; %1, 5, 7,10, 12 GPa
line_smbl={'-','--',':','-.'};
% figure('Position',[0.3 0.3 0.45 0.45],'Units','normalized');

% which grain growth laws to plot (choose max. 4)
% 1= plot, 0 = dont plot
K=1; %Karato
FS=0; %Faul and Scott
FJ=0; %Faul and Jackson
TH_fEn3=0; %Tasaka and Hiraga (3% En)
TH_fEn9=1; %Tasaka and Hiraga (9% En)
TH_fEn24=0; %Tasaka and Hiraga (24% En)
FF_fEn13_1GPa=1; % This study at 1GPa
FF_fEn13_7GPa=0; % This study at 7GPa
FF_fEn13_12GPa=1; % This study at 12GPa

% set a counter
nPlotCounter=0;

% plot grain growth equations for the different flow laws at different T
for i=1:length(T)
    T_K=T(i)+273.15; % Temperature(K)

    % (karato 1989/2012)
    if K
        nPlotCounter=nPlotCounter+1;
        n=2; %grain size exponent
        E=200*10^3;% Activation Energy(J/mol)
        k0=(10^4.6); % Pre exponential term (um^n s^-1)

        % Grain growth equation
        k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        %k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant (P dependent)

        GSt=((GS0_um^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt
        GSt_um = GSt;
%         GSt_um=GSt*1e6;

%         lgd=sprintf('%.0f ºC (K (1989)) ',T(i));
        lgd=sprintf('Karato (2008)');

        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
%         if i==length(T)
%             text(t(end),GSt_um(end),'Karato (2008)')
%         end

        hold on
    end
    
    %sol gel-Faul Scott 2006
    if FS
        nPlotCounter=nPlotCounter+1;
        n=4.3; %grain size exponent 
        E=390*10^3;% Activation Energy(J/mol)
        k0= 8.2e-14; % Pre exponential term (m^n s^-1)

        % Grain growth equation
        k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        %k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('Faul & Scott (2006)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;

%         if i==length(T)
%             text(t(end),GSt_um(end),'Faul & Scott (2006)')
%         end

        hold on
    end
    
    %sol gel-Faul Jackson 2007
    if FJ
        nPlotCounter=nPlotCounter+1;
        n=3.3; %grain size exponent 
        E=400*10^3;% Activation Energy(J/mol)
        k0= 1e-9; % Pre exponential term (m^n s^-1)

        % Grain growth equation
        k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        %k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('Faul & Jackson (2007)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
%         if i==length(T)
%             text(t(end),GSt_um(end),'Faul & Jackson (2007)')
%         end

        hold on
    end
    
    % Tasaka & Hiraga  (2013), Hiraga EPSL 2010
    %3 vol.%
    if TH_fEn3
        nPlotCounter=nPlotCounter+1;

        fEn=0.03; % Enstatite Fraction
        GS0_um=GS0*1e6;

        % Grain growth equation
        n=4; %grain size exponent
        
        p=polyfit([1360, 1310, 1260],[-2.02, -2.64, -3.45],1);% linear regression for (T,k)
        k =  p(1) * T(i) + p(2);%  (k at Temperature = T(i))
        k=10^k;
        
        GSt=((GS0_um^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt;%*1e6;

        lgd=sprintf('%.0f vol.%% En, Tasaka & Hiraga (2013) ',100*fEn);
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
        hold on
%         if i==length(T)
%             text(t(end),GSt_um(end),'(f_{En}=0.03), T&K(2013)','FontSize',18)
%         end
    end
    %9 vol.%
    if TH_fEn9
        nPlotCounter=nPlotCounter+1;

        fEn=0.09; % Enstatite Fraction
        GS0_um=GS0*1e6;
        % Grain growth equation
        n=4; %grain size exponent
        
        p=polyfit([1360, 1310, 1260],[-3.66, -4.28, -4.88],1);% linear regression for (T,k)
        k =  p(1) * T(i) + p(2); % x=Temperature in C
        k=10^k;
        
        GSt=((GS0_um^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt;%*1e6;

        lgd=sprintf('%.0f vol.%% En, Tasaka & Hiraga (2013) ',100*fEn);
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
        hold on
%         if i==length(T)
%             text(t(end),GSt_um(end),'(f_{En}=0.09), T&K(2013)','FontSize',18)
%         end
    end
    %24 vol.%
    if TH_fEn24
        nPlotCounter=nPlotCounter+1;

        fEn=0.24; % Enstatite Fraction
        GS0_um=GS0*1e6;
        % Grain growth equation
        n=4; %grain size exponent
        
        % Extrapolate k to T(i) at 24% En
        a=[];
        %3 %
        p=polyfit([1360, 1310, 1260],[-2.02, -2.64, -3.45],1);% linear regression for (T,k)
        a=[a,p(1)];
        %9 %
        p=polyfit([1360, 1310, 1260],[-3.66, -4.28, -4.88],1);% linear regression for (T,k)
        a=[a,p(1)];

        % 24 %
        % [1360],[-4.18]
        p = polyfit([3,9], a, 1);% linear regression
        a_24 =  p(1) * 24+ p(2);% interp a
        b_24 = -4.18-(a_24*1360);% get intercept at [1360],[-4.18]
        k = a_24*T(i)+b_24; % extrapolate to T(i) ;
        k=10^k;

        GSt=((GS0_um^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt;%*1e6;

        lgd=sprintf('%.0f vol.%% En, Tasaka & Hiraga (2013) ',100*fEn);
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
        hold on
%         if i==length(T)
%             text(t(end),GSt_um(end),'(f_{En}=0.09), T&K(2013)','FontSize',18)
%         end
    end
    
    % This study (1 GPa, 1400oC)
    if FF_fEn13_1GPa
        nPlotCounter=nPlotCounter+1;
        n=3.88; %grain size exponent 
        E=608*10^3;% Activation Energy(J/mol)
        k0= 2.11e-7; % Pre exponential term (m^n s^-1)
        V=4.3e-6; % Activation volume (m3/mol)
        P = [1]*1e9; %Pressure in Pa
        % Grain growth equation
%         k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('This study (1 GPa)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;

%         if i==length(T)
%             text(t(end),GSt_um(end),'This study')
%         end

        hold on
    end
     
     % This study (7 GPa, 1400oC)
    if FF_fEn13_7GPa
        nPlotCounter=nPlotCounter+1;
        n=3.88; %grain size exponent 
        E=608*10^3;% Activation Energy(J/mol)
        k0= 2.11e-7; % Pre exponential term (m^n s^-1)
        V=4.3e-6; % Activation volume (m3/mol)
        P = [7]*1e9; %Pressure in Pa
        % Grain growth equation
%         k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('This study (1 GPa)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;

%         if i==length(T)
%             text(t(end),GSt_um(end),'This study')
%         end

        hold on
     end
     
    % This study (12 GPa, 1400oC)
    if FF_fEn13_12GPa
        nPlotCounter=nPlotCounter+1;
        n=3.88; %grain size exponent 
        E=608*10^3;% Activation Energy(J/mol)
        k0= 2.11e-7; % Pre exponential term (m^n s^-1)
        V=4.3e-6; % Activation volume (m3/mol)
        P = [12]*1e9; %Pressure in Pa
        % Grain growth equation
%         k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (m)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('This study (12 GPa)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;

%         if i==length(T)
%             text(t(end),GSt_um(end),'This study')
%         end

        hold on
     end
    
end

%
%Plot experimental data
%

gg_data=load('C:\Users\filip\Documents\DOC\GBCD_Olivine\paper1\ggEXP\grainGrowthData.mat');
headers=gg_data.headers;
sampleNames=gg_data.sampleNames;
gg_data=gg_data.data;
gg_data=gg_data(:, 1:end-3);
P_exp=gg_data(1,:);
T_exp=gg_data(2,:);
t_exp=gg_data(3,:);
modeFit=gg_data(6,:);
meanfEn=gg_data(9,:);
nGrains=gg_data(10,:);

% Scatter plot data
hold on
for jj=1:size(gg_data,2)
    if P_exp(jj)>=1 && T_exp(jj)==1400
        if meanfEn(jj)<0.07
            scatter(t_exp(jj), modeFit(jj),100,mk_smbl{P_exp(jj)==P_all},'MarkerEdgeColor',cmapBlue2Red(find(T_exp(jj)==T),:),'MarkerEdgeAlpha',.6,'HandleVisibility','off')
        else 
            scatter(t_exp(jj), modeFit(jj),100,mk_smbl{P_exp(jj)==P_all},'filled','MarkerFaceColor',cmapBlue2Red(find(T_exp(jj)==T),:),'MarkerEdgeColor',cmapBlue2Red(find(T_exp(jj)==T),:),'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6,'HandleVisibility','off')
        end
    end
    hold on
end
% make symbol Legend
for jj=1:length(P_all)
    scatter(nan, nan, mk_smbl{jj},'MarkerEdgeColor','k','MarkerEdgeAlpha',.6,'DisplayName', [num2str(P_all(jj)), ' GPa'])
end
    scatter(nan, nan, mk_smbl{1},'MarkerEdgeColor','k','MarkerEdgeAlpha',.6,'DisplayName', '6 vol. % Px')
    scatter(nan, nan, mk_smbl{1},'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerEdgeAlpha',.6,'DisplayName', '13 vol. % Px')

hold off
legend('Location','bestOutside')

xlabel('Time [hours]','FontSize',16)
ylabel('Grain size [μm]','FontSize',16)
set(gca,'FontSize',14)
axis square
% f=gcf;f.Units='normalized';f.Position=[0.3 0.3 0.45 0.45]; 
set(gca,'yscale','log');
set(gca,'xscale','log'); % scale the x-axis

text(0.01,0.95, '1400 ºC','Units','normalized','FontSize',16)
%
%************-------------**************
%                      Sub plot 2                              
%************-------------**************

% Plot extrapolated grain growth laws to 10^n years:
subplot (2,1,2)
% figure('Units','normalized','Position', [0.05 0.05 0.8 0.5])
% Constants
% T=[1050,1200,1400,1520]; % Temperature (ºC)
T=[1400];
    
t_min=-1;% minimum time (10^n years)
t_max=9;% maximum time (10^n years)
t_resolution= 50;% number of points for each line
t=logspace(t_min,t_max, t_resolution); % %time points (years)
t_s=t*365*24*3600; % Time(s)

GS0=2.25*1e-6 ; % Initial grain size (m) 

lgd={};

% set a line style counter
nPlotCounter=0;

% plot grain growth equations for the different flow laws at different T
for i=1:length(T)
    T_K=T(i)+273.15; % Temperature(K)

    % (karato 1989/2012)
    if K
        nPlotCounter=nPlotCounter+1;
        n=2; %grain size exponent
        E=200*10^3;% Activation Energy(J/mol)
        k0=(10^4.6); % Pre exponential term (um^n s^-1)

        % Grain growth equation
        k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        %k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant (P dependent)

        GSt=((GS0_um^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
%         GSt_um=GSt*1e6;
        GSt_um = GSt;

%         lgd=sprintf('%.0f ºC (K (1989)) ',T(i));
        lgd=sprintf('Karato (2008)');

        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
%         if i==length(T)
%             text(t(end),GSt_um(end),'Karato (2008)')
%         end

        hold on
    end
    %sol gel-Faul Scott 2006
    if FS
        nPlotCounter=nPlotCounter+1;
        n=4.3; %grain size exponent 
        E=390*10^3;% Activation Energy(J/mol)
        k0= 8.2e-14; % Pre exponential term (m^n s^-1)

        % Grain growth equation
        k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        %k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('Faul & Scott (2006)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;

%         if i==length(T)
%             text(t(end),GSt_um(end),'Faul & Scott (2006)')
%         end

        hold on
    end  
    %sol gel-Faul Jackson 2007
    if FJ
        nPlotCounter=nPlotCounter+1;
        n=3.3; %grain size exponent 
        E=400*10^3;% Activation Energy(J/mol)
        k0= 1e-9; % Pre exponential term (m^n s^-1)

        % Grain growth equation
        k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        %k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('Faul & Jackson (2007)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
%         if i==length(T)
%             text(t(end),GSt_um(end),'Faul & Jackson (2007)')
%         end

        hold on
    end
    
    % Tasaka & Hiraga  (2013), Hiraga EPSL 2010
    %3 vol.%
    if TH_fEn3
        nPlotCounter=nPlotCounter+1;

        fEn=0.03; % Enstatite Fraction
        GS0_um=GS0*1e6;

        % Grain growth equation
        n=4; %grain size exponent
        
        p=polyfit([1360, 1310, 1260],[-2.02, -2.64, -3.45],1);% linear regression for (T,k)
        k =  p(1) * T(i) + p(2);%  (k at Temperature = T(i))
        k=10^k;
        
        GSt=((GS0_um^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt;%*1e6;

        lgd=sprintf('%.0f vol.%% En, Tasaka & Hiraga (2013) ',100*fEn);
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
        hold on
%         if i==length(T)
%             text(t(end),GSt_um(end),'(f_{En}=0.03), T&K(2013)','FontSize',18)
%         end
    end
    %9 vol.%
    if TH_fEn9
        nPlotCounter=nPlotCounter+1;

        fEn=0.09; % Enstatite Fraction
        GS0_um=GS0*1e6;
        % Grain growth equation
        n=4; %grain size exponent
        
        p=polyfit([1360, 1310, 1260],[-3.66, -4.28, -4.88],1);% linear regression for (T,k)
        k =  p(1) * T(i) + p(2);% x=Temperature
        k=10^k;
        
        GSt=((GS0_um^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt;%*1e6;

        lgd=sprintf('%.0f vol.%% En, Tasaka & Hiraga (2013) ',100*fEn);
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
        hold on
%         if i==length(T)
%             text(t(end),GSt_um(end),'(f_{En}=0.09), T&K(2013)','FontSize',18)
%         end
    end
    %24 vol.%
    if TH_fEn24
        nPlotCounter=nPlotCounter+1;

        fEn=0.24; % Enstatite Fraction
        GS0_um=GS0*1e6;
        % Grain growth equation
        n=4; %grain size exponent
        
        % Extrapolate k to T(i) at 24% En
        a=[];
        %3 %
        p=polyfit([1360, 1310, 1260],[-2.02, -2.64, -3.45],1);% linear regression for (T,k)
        a=[a,p(1)];
        %9 %
        p=polyfit([1360, 1310, 1260],[-3.66, -4.28, -4.88],1);% linear regression for (T,k)
        a=[a,p(1)];

        % 24 %
        % [1360],[-4.18]
        p = polyfit([3,9], a, 1);% linear regression
        a_24 =  p(1) * 24+ p(2);% interp a
        b_24 = -4.18-(a_24*1360);% get intercept at [1360],[-4.18]
        k = a_24*T(i)+b_24; % extrapolate to T(i) ;
        k=10^k;

        GSt=((GS0_um^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt;%*1e6;

        lgd=sprintf('%.0f vol.%% En, Tasaka & Hiraga (2013) ',100*fEn);
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;
        hold on
%         if i==length(T)
%             text(t(end),GSt_um(end),'(f_{En}=0.09), T&K(2013)','FontSize',18)
%         end
    end
    
   % This study (1 GPa, 1400oC)
    if FF_fEn13_1GPa
        nPlotCounter=nPlotCounter+1;
        n=3.88; %grain size exponent 
        E=608*10^3;% Activation Energy(J/mol)
        k0= 2.11e-7; % Pre exponential term (m^n s^-1)
        V=4.3e-6; % Activation volume (m3/mol)
        P = [1]*1e9; %Pressure in Pa
        % Grain growth equation
%         k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('This study (1 GPa)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;

%         if i==length(T)
%             text(t(end),GSt_um(end),'This study')
%         end

        hold on
    end
     
     % This study (7 GPa, 1400oC)
    if FF_fEn13_7GPa
        nPlotCounter=nPlotCounter+1;
        n=3.88; %grain size exponent 
        E=608*10^3;% Activation Energy(J/mol)
        k0= 2.11e-7; % Pre exponential term (m^n s^-1)
        V=4.3e-6; % Activation volume (m3/mol)
        P = [7]*1e9; %Pressure in Pa
        % Grain growth equation
%         k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('This study (1 GPa)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;

%         if i==length(T)
%             text(t(end),GSt_um(end),'This study')
%         end

        hold on
     end
      
    % This study (12 GPa, 1400oC)
    if FF_fEn13_12GPa
        nPlotCounter=nPlotCounter+1;
        n=3.88; %grain size exponent 
        E=608*10^3;% Activation Energy(J/mol)
        k0= 2.11e-7; % Pre exponential term (m^n s^-1)
        V=4.3e-6; % Activation volume (m3/mol)
        P = [12]*1e9; %Pressure in Pa
        % Grain growth equation
%         k=k0*exp(-(E)/(R*T_K)); % Rate Constant
        k=k0*exp(-(E+P*V)/(R*T_K)); % Rate Constant  (P dependent)

        GSt=((GS0^n)+(k.*t_s)).^(1/n); % Final Grain Size (um)
        GSt_um=GSt*1e6;% m to um

        lgd=sprintf('This study (12 GPa)');
        plot(t,GSt_um,line_smbl{nPlotCounter},'linewidth',2,'DisplayName',lgd,'Color', cmapBlue2Red_alpha(i,:))
        ax = gca; ax.ColorOrderIndex = i;

%         if i==length(T)
%             text(t(end),GSt_um(end),'This study')
%         end

        hold on
    end
     
end

hold off
legend('Location','bestOutside')

xlabel('Time [years]','FontSize',16)
ylabel('Grain size [μm]','FontSize',16)
set(gca,'FontSize',14)
axis square
% f=gcf;f.Units='normalized';f.Position=[0.3 0.3 0.45 0.45]; 
set(gca,'yscale','log');
set(gca,'xscale','log'); % scale the x-axis

text(0.01,0.95, '1400 ºC','Units','normalized','FontSize',16)

x_val=logspace(0,9,4);
        xticks(x_val)
        xticklabels(sprintfc('10^{%d}',log10(x_val)))

% lettersIds=fliplr({'a)','b)'});
% axs = findobj(f, 'type', 'axes', '-or', 'type', 'matlab.graphics.axis.Axes');
% 
% for axxN=1:length(lettersIds)
%     box(axs(axxN),'on');
%     t=title(axs(axxN),lettersIds(axxN),'FontSize',13);
%     set(t, 'horizontalAlignment', 'left')
%     set(t, 'units', 'normalized')
%     t.Position(1)=0;
% end
axis tight
        
print('C:\Users\filip\Documents\DOC\GBCD_Olivine\paper1\figures\DataExtrapolation_new','-dpdf', '-painters')
% saveFigure('/Users/fereira/Documents/DOC/GBCD_Olivine/paper1/figures/DataExtrapolation.png')



%%
function makeScatterLegend
    hold on
    s_md=scatter(nan, nan,'ks', 'DisplayName',' Mode');
    hold on
    s_md=scatter(nan, nan,'kd', 'DisplayName',' Median');
    hold on
    s_m=scatter(nan, nan,'ko','DisplayName',' Mean');
    hold off
end

function cmapBlue2Red=cmap(length_cmap)
    blue = [0 0 1]; %[29, 173, 207]/255;
    red = [1 0 0]; %[249, 84, 28]/255;
    cmapBlue2Red = [linspace(blue(1),red(1),length_cmap)', linspace(blue(2),red(2),length_cmap)', linspace(blue(3),red(3),length_cmap)'];
end

function setLogXticks
        x_values=[0 1 2 5 10 20 30 40 70 110];
        xticks(x_values)
        xticklabels(sprintfc('%d',x_values))
end

function [ebsd_clean,grains_clean]=ebsdCleaner(ebsd)
% Cleaning parameters
    % Minimum confidence index
    minCI =0.01; 
    % Minimum  grain size (points)
    minGS =20; 
    % Minimum grain misorientation angle threshold
    minANG =10*degree; 
    % Grain boundaries smooth factor (number of iterations)
    sF=10;
    % ------Clean data-----------
        % Remove data with CI lower than the defined minimum
       
        idNotIdx=0;
        

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
         % Smooth grain boundaries
        grains_clean=smooth(grains,sF);
        ebsd_clean=ebsd;
end
