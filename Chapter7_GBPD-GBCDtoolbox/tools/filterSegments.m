function n=filterSegments(segments_fname_in, segProperty,segValue,segments_fname_out)
%Function actions:
    %Filter segments that have a given property (segProperty) within a given
    %value range(segProperty)
%Input:
    %segments_fname_in = Path to the segments input file (including extension)
    % segProperty = Property of interest (e.g. 'misAngle', 'segLength','traceAngle', 'numberOfSegments', TODO: misAxis, orientation)
    % segValue = array with min and max values e.g. [min max]. if nan([nannan]), calc min/max
    % segments_fname_out=
    % GBCD_par = gbpd/gbcd/parameters (obtained from getGBCDpar.m)
%Output:
    % File created (segments_fname_out) with the filtered segments
    % n = Number of filtered segments
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
 
  % Read all segments
    A = regexp(fileread(segments_fname_in), '\n', 'split'); % read all segments
    A=cellfun(@strtrim, A,'UniformOutput' ,false); %remove empty spaces
    A=A(~cellfun('isempty',A)); % remove empty lines
    
    %Number of header lines
        hL=0;
        for nline=1:50 %check first 50 lines
            if ischar(A{nline}) && startsWith(A{nline},'#')
              hL = hL+1;
            end
        end
        
        fprintf('\n ')
      
 if any(size(segValue)==2)
     
     segm = readtable(segments_fname_in,'CommentStyle', '#' );
    switch segProperty
        case 'misAngle'
            prop=segm.Var7; % misorientation angle  in degrees
        case 'segLength'
            prop=segm.Var14; % segment length in um
        case 'traceAngle'
            prop=segm.Var15; % trace angle  in degrees
        otherwise
            error([segProperty, ' is not a known segments property.'])
    end
        clear segm; % huge file so clear memory after use

        % if value is nan use min or max
        if isnan(segValue(1))
            segValue(1)=min(prop);
        end
        if isnan(segValue(2))
            segValue(2)=max(prop);
        end

        % get id of filtered segments
        id=(prop>segValue(1) & prop<=segValue(2));
        n= sum(id);
        disp(['Filtered ' ,num2str(n),' out of ' ,num2str(length(id)), ' segments (' ,num2str(round(100*n/length(id))),'%)'])
        id =[true(hL,1); id]; % keep header lines
        A=A(id); % remove not wanted lines
 elseif size(segValue)==1% single value
     switch segProperty
        case 'numberOfSegments'
            s = RandStream('mlfg6331_64','seed','shuffle');
            if segValue<(length(A)-hL)
                id=datasample(s, (hL+1):length(A), segValue, 'Replace', false)';
            else
                id=datasample(s, (hL+1):length(A), segValue, 'Replace', true)';
                str=sprintf('More values to sample (%.0f) than actual data (%.0f). Some segments will be repeated.',segValue,(length(A)-hL));
                warning(str)
            end
            n= length(id);
            disp(['Filtered ' ,num2str(n),' out of ' ,num2str(length(A)-hL), ' segments (' ,num2str(round( 100*n/(length(A)-hL) )),'%)'])
            A=[{A{1:hL}}';{A{id}}'];
         otherwise
            error([segProperty, ' is not a known segments property.'])
     end
 end
     
    % Write Filtered segments
    fid = fopen(segments_fname_out, 'w'); % open file in write mode
    fprintf(fid, '%s\n', A{:}); % print all segments to file
    fclose(fid);
    
fprintf('\n ')