function CS=getCS(ebsdFname)
%Function actions:
    % Get crystal symmetries from .ang file
%Input:
    % ebsdFname = Path to EBSD file
%Output:
    % CS = Cell with crystal symmetries
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

% read file header
    
    % colors for EBSD phase plotting
    EBSDColorNames = {'light blue','light green','light red',...
    'cyan','magenta','yellow',...
    'blue','green','red',...
    'dark blue','dark green','dark red'};
  
hl = file2cell(ebsdFname,2000);
phasePos = strmatch('# Phase ',hl);

  if isempty(phasePos)
    phasePos = strmatch('# MaterialName ',hl)-1;
  end
        
    for i = 1:length(phasePos)
      pos = phasePos(i);
      
      % load phase number
      phase = readByToken(hl(pos:pos+100),'# Phase',i);
            
      % load mineral data
      mineral = readByToken(hl(pos:pos+100),'# MaterialName');
      laue = readByToken(hl(pos:pos+100),'# Symmetry');
      lattice = readByToken(hl(pos:pos+100),'# LatticeConstants',[1 1 1 90 90 90]);
      
      % setup crytsal symmetry
      options = {};
      switch laue
        case {'-3m' '32' '3' '62' '6'}
          options = {'X||a'};
        case '2'
          options = {'X||a*'};
        case '1'
          options = {'X||a'};
        case '20'
          laue = {'2'};
          options = {'y||b','z||c'};
        otherwise
          if lattice(6) ~= 90
            options = {'X||a'};
          end
      end
      CS{phase+1} = crystalSymmetry(laue,lattice(1:3)',lattice(4:6)'*degree,'mineral',mineral, 'color',EBSDColorNames{phase}, options{:});
    end
    CS{1}='notIndexed' ;
function value = readByToken(cellStr,token,default)

    values = regexp(cellStr,[token '\s*(.*)'],'tokens');
    id = find(~cellfun(@isempty,values),1);
    if ~isempty(id)
        value = strtrim(char(values{id}{1}));

        if nargin > 2 && ~isempty(default) && isnumeric(default)
            value = str2num(value);
        end

    elseif nargin > 2
        value = default;
    else 
        value = [];
    end

end
end