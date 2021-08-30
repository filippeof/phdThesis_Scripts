function  plotGBCD(gbcds, cs, axs, angs,varargin)
%Function actions:
    % Plot GBCD in a stereographic projection
%Input:
    % gbcds = GBCD structure(s) containing rho, theta and m.u.d.
    % ax= Axis of misorientation (MTEX class Miller)
    % ang= Angles of misorientation in degrees (double array)
%Output:
    % Stereographic projection plot of the GBCD
%Flags:
    % 'scatter' : Plot scatter plot of the GBCD
    % 'smooth' : Interpolated plot
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
%% Plot gbcd
fprintf(1,'\n -plotGBCD: Plotting GBCD...')

mtexFig=newMtexFigure;


for i=1:length(gbcds)
    gbcd=gbcds{i};
    ax=round(Miller(axs(i),'uvw'));
    ang=angs(i);
   
    if i~=1
        mtexFig.nextAxis;
    else
        mtexFig.gca;
    end
    
    vec=vector3d('polar', deg2rad(gbcd.rho),deg2rad(gbcd.theta));
    
% Scatter plot   
    if contains(varargin, 'scatter')
        if contains(varargin, 'fundamentalRegion')
            m=Miller(vec,cs);
            plot(m,gbcd.mud,'fundamentalRegion','parent', mtexFig.gca)
        else
            plot(vec,gbcd.mud,'upper','parent', mtexFig.gca)
        end
% Interpolation plot   
    else 
        if contains(varargin, 'fundamentalRegion')
            m=Miller(vec,cs);
            plot(m,gbcd.mud,'fundamentalRegion','smooth','parent', mtexFig.gca)
        else
            grid_res=round(mode(diff(gbcd.theta))); % Grid resolution: angle interval between values
            equiGrid=plotS2Grid('resolution',grid_res*degree);

            idwVec = interp(vec,gbcd.mud,equiGrid,'inverseDistance'); 
            pcolor(equiGrid,idwVec,'projection','edist','upper','smooth','parent', mtexFig.gca); % 
        end
    
    end
    % Set title    
    strTitle= sprintf('[%.0f %.0f %.0f] %.0f%c', round(ax.u),round(ax.v),round(ax.w),ang,char(176));
    mtexTitle(strTitle,'interpreter','tex','parent', mtexFig.gca)%

end

fprintf(1,' Done! \n ')
end