function  plotGBPD(gbpd,cs,varargin)
%Function actions:
    % Plot GBPD in a stereographic projection
%Input:
    % gbpd = GBPD structure containing rho, theta and m.u.d.
%Output:
    % Stereographic projection plot of the GBPD
%Flags:
    % 'scatter' : Plot scatter plot of the GBPD
    % 'fundamentalRegion' : Plot in the fundamental Region of the crystal symmetry 
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

%% Plot gbpd
fprintf(1,'\n -plotGBPD: Plotting GBPD...')

vec=vector3d('polar', deg2rad(gbpd.rho), deg2rad(gbpd.theta));

% Scatter plot   
if ~isempty(varargin) && any(contains(varargin, 'scatter'))
      if contains(varargin, 'fundamentalRegion')
        m=Miller(vec,cs);
        figure;plot(m,gbpd.mud,'fundamentalRegion')
    
    else
         figure;plot(vec,gbpd.mud,'upper')
      end
% Interpolation plot   
else 
    if ~isempty(varargin) && any(contains(varargin, 'fundamentalRegion'))
        m=Miller(vec,cs);
        figure;plot(m,gbpd.mud,'fundamentalRegion','smooth')
    else
        grid_res=round(mode(diff(gbpd.theta))); % Grid resolution: angle interval between values (in deg)
        % grid_res=9;
        equiGrid=plotS2Grid('resolution',grid_res*degree);

        idwVec = interp(vec,gbpd.mud,equiGrid,'inverseDistance'); 
        pcolor(equiGrid,idwVec,'projection','edist','upper','smooth'); % 
    end

end
% Plot Miller index above it
% h = [cs.aAxis,cs.bAxis,cs.cAxis];
% hold on
% plot(h,'markerFaceColor',[ 0.5765 0.2314 0.2549],'markerSize',5, 'labeled','all', 'fontSize',12,'Color',[ 0.5765 0.2314 0.2549],'upper')
% hold off
%
% set(gca,'fontSize',18)
% 
% mtexColorMap(parula)

% hold on
%     plotMirrorPlanes(cs,'upper')
% hold off

%% Plot max/min

% [~,id_max_mud]=max(gbpd.mud);
% [~,id_min_mud]=min(gbpd.mud);
% 
% % m_max=round(Miller(vec(id_max_mud),cs));
% % m_min=round(Miller(vec(id_min_mud),cs));
% m_max=Miller(vec(id_max_mud),cs);
% m_min=Miller(vec(id_min_mud),cs);
% 
%%Low id Miller idx
% maxIdx=4;
% equiGrid=plotS2Grid('resolution', 2*degree,'fundmentalRegion'); 
% m=Miller(equiGrid,cs);
% m=round(m); m=unique(m);
% m=m(m.h<=maxIdx & m.k<=maxIdx & m.l<=maxIdx);
% m=symmetrise(m,cs,'upper');m=m(:);
% 
% [~,id_max]=min(angle(m_max,m));
% m_max=m(id_max);
% [~,id_min]=min(angle(m_min,m));
% m_min=m(id_min);
% % f=newMtexFigure;plot(m, 'labeled','all', 'fontSize',8,'upper','figSize','huge','markerSize',3,'grid','grid_res',pi/2)
% 
% %plot Miller IDX
% hold on
% plot(m_max,'markerFaceColor',[0 0 0],'markerSize',8, 'labeled','all', 'fontSize',16,'Color',[0 0 0],'upper','parent',gca)
% hold off
% 
% hold on
% plot(m_min,'markerFaceColor',[1 1 1],'markerSize',8, 'labeled','all', 'fontSize',16,'Color',[1 1 1],'upper','parent',gca)
% hold off

% Correct colorbar position
% set(gca,'fontSize',18)
% 
% mtexColorMap(parula)
% vbar=mtexColorbar('fontSize', 18);
% ylabel(vbar,'m.u.d.','Interpreter','tex','fontSize',18);
% vbarPos = get(vbar,'Position');
% x_offset=round(vbarPos(1)*0.05);
% set(vbar,'Position',[vbarPos(1)+x_offset vbarPos(2) vbarPos(3) vbarPos(4)])
fprintf(1,' Done! \n ')
end