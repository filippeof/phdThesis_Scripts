function n=numberOfSegments(segmentsFname)
%Function actions:
    % Get number of segments in a file
%Input:
    % segmentsFname = Path to the segments file
%Output:
    % n = Number of segments

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
        A = regexp(fileread(segmentsFname), '\n', 'split')';
        A=cellfun(@strtrim, A,'UniformOutput' ,false); %remove empty lines
        A=A(~cellfun('isempty',A)); %remove empty lines

        %Number of header lines
        hL=0;
        for nline=1:50 %check first 50 lines
            if ischar(A{nline}) && startsWith(A{nline},'#')
              hL = hL+1;
            end
        end
        
        n=length(A)-hL;