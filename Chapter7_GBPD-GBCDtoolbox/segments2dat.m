function datFname = segments2dat(GBCDfolder,segments_fname,GBCD_par)
%Function actions:
    %1- Copy boundary segments file to GBPD folder 
    %2- Change  input.txt file
    %3- run GBPD program
%Input:
    % GBCDfolder = root folder with gbpd/gbcd programs
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
%% Set file names, folder paths
    oldFolder=pwd;
    cd (GBCDfolder)
    [~,sample_name,ext] = fileparts(segments_fname) ;
    % output file name should be:
     datFname = fullfile(GBCDfolder, 'graph_pd',['gbpd_', sample_name, '_2d_gmt.dat']);

 %% Replace location of boundary segments file in 'calc GBPD' input file
    calcGbpdFolder = fullfile(GBCDfolder,'calc_gbpd');
    inputFname=fullfile(calcGbpdFolder,'input.txt');

    % replace second line of input file with the name of segment file
    clear A; A = regexp(fileread(inputFname), '\n', 'split');% read all data to A, line by line
    A{2} = sprintf([sample_name,ext]);% replace content of 2nd line

    % Replace parameters in 'calc GBPD' input file
    %GBCD_par=[msym, rot, rad, resolution, seg, stepSize, cutOff];
    A{4} = sprintf( '%d\t%d\t%d', GBCD_par(1:3)); % replace symmetry, rotation, radian
    A{6} = sprintf( '%d\t%d', GBCD_par(4:5)); % replace resolution
    A{8} = sprintf( '%d', GBCD_par(6)); % replace segment analysis
    A{10} = sprintf( '%.2f', GBCD_par(7)); % replace step size
    A{12} = sprintf( '%.2f', GBCD_par(8)); % replace small length cutoff

    fid = fopen(inputFname, 'w'); % Open 'input.txt'
    fprintf(fid, '%s\n', A{:}); %
    fclose(fid);

    %copy file
    copyfile(segments_fname,fullfile(calcGbpdFolder,[sample_name,ext]) , 'f');

%% Replace parameters in  'plot GBPD' input file
    inputFname=fullfile(GBCDfolder,'graph_pd','input.txt');

    clear A; A = regexp(fileread(inputFname), '\n', 'split');% read all data to A, line by line

    % replace file name 
    A{1} = sprintf(['gbpd_',sample_name,ext]);% replace content of 1st line

    %GBCD_par=[msym, rot, rad, resolution, seg, stepSize, cutOff];
    A{3} = sprintf( '%d ', GBCD_par(1));% replace symmetry
    A{5} = sprintf( '%d\t%d', GBCD_par(4:5));% replace resolution

    fid = fopen(inputFname, 'w');% open input.txt
    fprintf(fid, '%s\n', A{:}); % write A 
    fclose(fid); % close

    %% Run GBPD script

    % Set Fortran path For MATLAB (if not already):
    try
        if contains(computer , 'MAC' ) % MAC
            if ~contains(getenv('PATH'), ':/usr/local/bin') 
                setenv('PATH', [getenv('PATH') ':/usr/local/bin']); %set path for mac
            end
            system(['chmod +x ',fullfile(GBCDfolder, 'gbpd.sh')]);% make shell file executable (Mac/Linux)
            system(['. ',fullfile(GBCDfolder, 'gbpd.sh')]);% Run and Show output 
            fprintf('\n -segments2dat: GBPD wrote to file: %s \n ', datFname)

        elseif contains(computer , 'LNX' ) % Linux
                system(['chmod +x ',fullfile(GBCDfolder, 'gbpd.sh')]);% make shell file executable (Mac/Linux)
                system(['. ',fullfile(GBCDfolder, 'gbpd.sh')]);% Run and Show output
                fprintf('\n -segments2dat: GBPD wrote to file: %s \n ', datFname)


        elseif contains(computer , 'WIN' )
            system(['chmod +x ',fullfile(GBCDfolder, 'gbpd_WIN.sh')]);% make shell file executable (Mac/Linux)
            system(['sh ',fullfile(GBCDfolder, 'gbpd_WIN.sh')]);% Run and Show output
            fprintf('\n -segments2dat: GBPD wrote to file: %s \n ', datFname)

        end
        
     catch
            error('It was not possible to execute the shell Files. ')
    end
   % [x,y]=system([GBCDfolder, '/gbpd.sh']) % work with output
    % [~, ~]=system([GBCDfolder, '/gbpd.sh']);% omit output
    %add GMT folder to path
    %setenv('PATH', [getenv('PATH') ':/opt/local/bin']); %set GMT path for mac
    %system('./Draw_stereograms 1 gbpd_SrTiO3_all_EDAX_2d_gmt 2d rainbow 0.76 1.37 0.06 stereo CUB')

    cd (oldFolder)
end
