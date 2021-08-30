function joinSegmentsNormalized(fnamesIn,fnameOut,varargin)
%Function actions:
    % Join multiple segment files into one
%Input:
    % fnameOut = Path to the output file (including extension)
    % fnamesIn =  Cell with paths to the input files (including extension)
    %
%Output:
    % Segments files (fnamesIn) joined together into a single file(fnameOut)
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
% Create random stream
    s = RandStream('mlfg6331_64');
% open output file in write mode(delete all content)
    fid = fopen(fnameOut, 'w'); 
    fclose(fid);
    
    %Normalize by user defined number of segments to extract
    if nargin==4 && contains(varargin, 'min')
        minLines=str2double(varargin{end});
    else
        % Normalize by minimum of segments among all files
        numlines=zeros(1,length(fnamesIn));
        % Reads the number of lines of each file
       for j=1:length(fnamesIn)
            A = regexp(fileread(fnamesIn{j}), '\n', 'split')';
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
    % open output file in append mode
    fid = fopen(fnameOut,'a'); 
    
    for j=1:length(fnamesIn)
        clear A A_new
         A = regexp(fileread(fnamesIn{j}), '\n', 'split')';% read all segments
         %Number of header lines
        hL=0;
        for nline=1:50 %check first 50 lines
            if ischar(A{nline}) && startsWith(A{nline},'#')
              hL = hL+1;
            end
        end
        
         if j==1 %if its the first file, keep headers
            A=cellfun(@strtrim, A,'UniformOutput' ,false); %remove empty lines
            A=A(~cellfun('isempty',A)); %remove empty lines
            idsA=datasample(s,(hL+1):length(A),minLines,'Replace',false)';
            A_new=[{A{1:hL}}';{A{idsA}}'];
         else
            A(1:hL)=[];% Remove headers of subsequent files
            A=cellfun(@strtrim, A,'UniformOutput' ,false); %remove empty lines
            A=A(~cellfun('isempty',A)); %remove empty lines
            idsA=datasample(s,1:length(A),minLines,'Replace',false)';
            A_new={A{idsA}};
         end
        fprintf(fid, '%s\n', A_new{:}); %print all segments to file
    end

    fclose(fid);
    fprintf('\n -joinSegmentsNormalize: %.0f segments were extracted from each of the %.0f files and joined to: %s \n ', minLines, length(fnamesIn), fnameOut)

    end
