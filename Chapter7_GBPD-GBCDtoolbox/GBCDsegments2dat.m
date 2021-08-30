function datFname = GBCDsegments2dat(ax,ang,GBCDfolder,segments_fname, GBCD_par)
%Function actions:
    %1- Copy boundary segments file to GBCD folder 
    %2- Change  input.txt file
    %3- run GBCD program
%Input:
    % ax= Axis of misorientation (MTEX class Miller)
    % ang= Angles of misorientation in degrees (double array)
    % GBCDfolder = root folder with gbpd/gbcd programs
    % segments_fname = path to segments file
    % GBCD_par = gbpd/gbcd/parameters (obtained from getGBCDpar.m)
%Output:
    %datFname: Path to the GBCD .dat file located in the gbcd_graph_fd folder to be ploted
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
%% Set file names, folder paths
oldFolder=pwd;
cd (GBCDfolder)

if length(ax)~=length(ang)
       error ('Different axis/ angle size. They must be the same size.')
end
[~,sample_name,ext] = fileparts(segments_fname) ;

% output file names should be:
for i=1:length(ang)
 datFname{i} = fullfile(GBCDfolder,'gbcd_graph_fd',['gbcd_',sample_name, '_gmt_',num2str(i),'.dat']);
end

%% Replace location of boundary segments file in 'calc GBCD' input file
inputFname = fullfile(GBCDfolder,'calc_gbcd_stereo_fd','input.txt');

% replace second line of 'gbcd calc' input file with the name of segment file
clear A; A = regexp(fileread(inputFname), '\n', 'split');% read all data to A, line by line
A{2} = sprintf([sample_name,ext]);% replace content of 2nd line
fid = fopen(inputFname, 'w');% open input.txt

% Replace parameters in 'calc GBCD' input file
%GBCD_par=[msym, rot, rad, resolution, seg, stepSize, cutOff];
A{4} = sprintf( '%d\t%d\t%d', GBCD_par(1:3)); % replace symmetry, rotation, radian
A{6} = sprintf( '%d\t%d', GBCD_par(4:5)); % replace resolution
A{8} = sprintf( '%d', GBCD_par(6)); % replace segment analysis
A{10} = sprintf( '%4.1f', GBCD_par(7)); % replace step size
A{12} = sprintf( '%4.2f', GBCD_par(8)); % replace small length cutoff

fprintf(fid, '%s\n', A{:});
fclose(fid);

%copy file
copyfile(segments_fname,fullfile(GBCDfolder,'calc_gbcd_stereo_fd', [sample_name,ext]), 'f');

%% Replace file path, axis and angle in  'plot GBCD' input file
inputFname = fullfile(GBCDfolder,'gbcd_graph_fd','input.txt');

% replace 12+(1:length(ang)) lines of 'gbcd plot' input file with the ax(x,y,z)/ang(deg)
clear A; A = regexp(fileread(inputFname), '\n', 'split');% read all data to A, line by line
% replace file name 
    A{1} = sprintf(['gbcd_',sample_name,ext]);% replace content of 1st line

A(13:end) = [];% clear existing ax/angle
A(~cellfun('isempty',A)); % clear empty lines

for i=1:length(ang)
    A{12+i} = sprintf( '%4.2f \t %4.2f \t %4.2f \t %4.1f', ax(i).x, ax(i).y, ax(i).z, ang(i)');% replace content of 13th to (12+max(i)) lines
end
% Replace parameters in 'plot GBCD' input file
%GBCD_par=[msym, rot, rad, resolution, seg, stepSize, cutOff];
A{3} = sprintf( '%d ', GBCD_par(1));% replace symmetry
A{9} = sprintf( '%d\t%d', GBCD_par(4:5));% replace resolution
A{11} = sprintf( '%d ', length(ax));% replace number of plots

fid = fopen(inputFname, 'w');% open input.txt
fprintf(fid, '%s\n', A{:}); % write A 
fclose(fid); % close

  %% Run GBCD script

    % Set Fortran path For MATLAB (if not already):
try
    if contains(computer , 'MAC' ) % MAC
        if ~contains(getenv('PATH'), ':/usr/local/bin') 
            setenv('PATH', [getenv('PATH') ':/usr/local/bin']); %set path for mac
        end
        system(['chmod +x ',fullfile(GBCDfolder, 'gbcd.sh')]);% make shell file executable (Mac/Linux)
        system(['. ',fullfile(GBCDfolder, 'gbcd.sh')]);% Run and Show output 
        fprintf('\n-GBCDsegments2dat: GBCD wrote to file(s): \n ')
        disp(datFname);
        fprintf('\n')

    elseif contains(computer , 'LNX' ) % Linux
            system(['chmod +x ',fullfile(GBCDfolder, 'gbcd.sh')]);% make shell file executable (Mac/Linux)
            system(['. ',fullfile(GBCDfolder, 'gbcd.sh')]);% Run and Show output
            fprintf('\n-GBCDsegments2dat: GBCD wrote to file(s): \n ')
            disp(datFname);
            fprintf('\n')

    elseif contains(computer , 'WIN' )
        system(['chmod +x ',fullfile(GBCDfolder, 'gbcd_WIN.sh')]);% make shell file executable (Mac/Linux)
        system(['sh ',fullfile(GBCDfolder, 'gbcd_WIN.sh')]);% Run and Show output
        fprintf('\n-GBCDsegments2dat: GBCD wrote to file(s): \n ')
        disp(datFname);
        fprintf('\n')
    end

 catch
        error('It was not possible to execute the shell Files. ')
end
%add GMT folder to path
%     setenv('PATH', [getenv('PATH') ':/opt/local/bin']); %set path for mac

cd (oldFolder)
end