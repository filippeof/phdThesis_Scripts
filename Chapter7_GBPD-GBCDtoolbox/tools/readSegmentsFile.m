function dataFile=readSegmentsFile(segmentsFname)
%Function actions:
    % Read segments file into a table.  Acces data afterwards with e.g. 
    % 'mis_angle= dataFile.m_angle'
%Input:
    % segmentsFname = Path to the segments file
%Output:
    % dataFile = Table with the segments data

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

% fname= Segments file name
% fname='/Users/fereira/Documents/DOC/GBCD_Olivine/paper3/paper3-Test/map20170924164021352_crop_clean_segments_mtex.txt';
opt={'CommentStyle', '#' , 'FileType','text','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines', 0, 'ReadVariableNames', false,'TreatAsEmpty',' '};
    dataFile = readtable(segmentsFname,opt{:});
   if size (dataFile,2) == 13 || size(dataFile,2) == 15 || size(dataFile,2) == 22% empty space at the last column
       if isempty(dataFile.(dataFile.Properties.VariableNames{end})(1))
            dataFile.(dataFile.Properties.VariableNames{end})=[];
       end
   end
   
    if size (dataFile,2) == 21
        dataFile.Properties.VariableNames = {'o_right_rad_Phi1','o_right_rad_PHI','o_right_rad_Phi2', 'o_left_rad_Phi1','o_left_rad_PHI','o_left_rad_Phi2',...
            'm_angle','m_axis_right_u','m_axis_right_v','m_axis_right_w','m_axis_left_u','m_axis_left_v','m_axis_left_w',...
            'trace_length','trace_angle','X1','Y1','X2','Y2','Id_grain_R','Id_grain_L'};
    
    elseif size (dataFile,2) == 14
        dataFile.Properties.VariableNames = {'o_right_rad_Phi1','o_right_rad_PHI','o_right_rad_Phi2', 'o_left_rad_Phi1','o_left_rad_PHI','o_left_rad_Phi2',...
            'trace_length','trace_angle','X1','Y1','X2','Y2','Id_grain_R','Id_grain_L'};
        
    elseif size (dataFile,2) == 12
        dataFile.Properties.VariableNames = {'o_right_rad_Phi1','o_right_rad_PHI','o_right_rad_Phi2', 'o_left_rad_Phi1','o_left_rad_PHI','o_left_rad_Phi2',...
            'trace_length','trace_angle','X1','Y1','X2','Y2'};
        
    else
        error('Segment file unknown')
    end
         fprintf('\n -readSegmentsFile: Segments data imported from file: \n %s ',segmentsFname)
         fprintf('\n ')
end
       

% To export table back to segmentsFile:
%
% segmentsFnameOut='testOut.txt';
% [fid] = fopen(segmentsFnameOut, 'w');% if file exists it will be overwritten
% 
% fprintf(fid,' %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f   %3.0f %3.0f %3.0f    %3.0f %3.0f %3.0f  %8.2f %8.2f    %8.2f    %8.2f   %8.2f   %8.2f   %4.0f   %4.0f\n',...
%   table2array(dataFile));
% fclose(fid);
