function joinSegments(fnameOut,fnamesIn)
%Function actions:
    % Join multiple segment files into one
%Input:
    % fnameOut = Path to the output file (including extension)
    % fnamesIn =  Cell with paths to the input files (including extension)
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
%Check if fnamesOut is in fnamesIn and remove it from the list
outputInsideList=ismember(fnamesIn,fnameOut);
if any(outputInsideList)
    fnamesIn(outputInsideList)=[];
end
% open output file in write mode (delete all content)
    fid = fopen(fnameOut, 'w'); 
    fclose(fid);
% open output file in append mode
     fid = fopen(fnameOut,'a');
    for j=1:length(fnamesIn)
        clear A A_new
         A = regexp(fileread(fnamesIn{j}), '\n', 'split')';% read all segments
         %Number of header lines
        hL=0;
        for nline=1:20 %check first 50 lines
            if ischar(A{nline}) && startsWith(A{nline},'#')
              hL = hL+1;
            end
        end
        
         if j==1 %if its the first file:
            A=cellfun(@strtrim, A,'UniformOutput' ,false); %remove empty lines
            A=A(~cellfun('isempty',A)); %remove empty lines
         else
            A(1:hL)=[];% Remove headers of subsequent files
            A=cellfun(@strtrim, A,'UniformOutput' ,false); %remove empty spaces
            A=A(~cellfun('isempty',A)); %remove empty lines
         end
        fprintf(fid, '%s\n', A{:}); %print all segments to file
    end
        fclose(fid);
fprintf('\n -joinSegments: %.0f files joined to: %s \n ',length(fnamesIn), fnameOut)

end