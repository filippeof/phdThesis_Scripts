function [gb,mp,rbv,grains]= slipTrev(grains,phaseName,section)
% calculate m' factor and residual burger vector for A, B, C, D, E -types olivine slip system
setMTEXpref('xAxisDirection','east');setMTEXpref('zAxisDirection','intoPlane');
%grains should be of single phase
gb=grains.boundary;
% id_gbOut=setdiff(grains.id,grains(phaseName).id); %id of grains that are not olivine
% [rowOut,~]=find(gb.grainId(:,:)==id_gbOut'); 
% gb_in=true(length(gb), 1);
% gb_in(rowOut)=false;
% 
% gb=gb(gb_in);
gb=gb(phaseName,phaseName);

ori=grains.prop.meanRotation;

% Define Slip systems [n; b]
slipS(:,:,1)=[ 0 1 0  ; 1 0 0];% a type
slipS(:,:,2)=[ 0 1 0  ; 0 0 1];% b-type
slipS(:,:,3)=[ 1 0 0  ; 0 0 1];% c-type
slipS(:,:,4)=[ 0 1 1  ; 1 0 0];% d-type
slipS(:,:,5)=[ 0 0 1  ; 1 0 0];% e type
slipS(:,:,6)=[ 1 1 0  ; 0 0 1];% add new ss
slipS(:,:,7)=[ 1 3 0  ; 0 0 1];% add new ss
% slipS(:,:,6)=[ n1 n2 n3  ; b1 b2 b3];% add new ss

nSS=size(slipS,3); % number of slip systems

cs=crystalSymmetry('mmm', [4.762 10.225 5.994], 'mineral', 'olivine');
CRSS=ones(nSS,1);% Same CRSS
clear n b SS
for i =1:nSS
    n{i}=Miller(slipS(1,1,i),slipS(1,2,i),slipS(1,3,i),cs,'hkl');% n=slip plane normal
    b{i}=Miller(slipS(2,1,i),slipS(2,2,i),slipS(2,3,i),cs,'uvw');%b=slip direction
    % check if angle b to n is 90 degrees (i.e. that slip direction is in the slip plane)
    Angle_n_to_b = round(angle(n{i},b{i})./degree);
    if Angle_n_to_b ==90
        SS{i} = slipSystem(b{i},n{i},CRSS(i));
    end
end
n=[n{:}]; b=[b{:}];SS=[SS{:}];

% arrange SS in symmetrical groups
id_Sym=[]; clear symmSS
SS_sym= SS.symmetrise;
for i=1:nSS
    SSi_sym=SS(i).symmetrise;
    is_SSSym = arrayfun(@(x)find(SSi_sym==x,1), SS_sym, 'UniformOutput', false);
    symmSS{i} =find(~cellfun(@isempty,is_SSSym)); % each cell has the ids of the symmetrical SS in SS_SYM
    id_Sym=[id_Sym; repmat(i, length(symmSS{i}), 1)];% Array with id of parent SS
end

    inv_ori=inv(ori);
    P=0.3; %GPa
    sig=0.168; %GPa
    sigma=stressTensor([P 0 0; 0 P sig ; 0 sig P]); % torsion in tensor notation modif. from Tielke 2016( rot Y, -90º)

    switch section
        case 'xz'
            g_ctd=grains.centroid;
            x=g_ctd(:,1);
            y=g_ctd(:,2);
            theta= atan(y./x);
            rot = rotation('axis',yvector,'angle',theta);
            sigma = rotate(sigma, rot); 
        case 'yz'
            rot = rotation('axis',zvector,'angle',pi/2);
            sigma = rotate(sigma, rot); 
        case 'xy'
            % nothing to do
    end
         
    sigmaLocal=rotate(sigma, inv_ori);% sigma into specimen coordinate
    SF= SS_sym.SchmidFactor(sigmaLocal);% schimidt factor
    [~, idSSactive]= max (SF,[],2);% Maximum along rows -> get id of active slip system
%% Plot active slip system
% grains.prop.activeSSid=idSSactive;% create a property for the id of the active slip system
% figure
grains.prop.activeSSid_sym=changem(idSSactive,id_Sym,1:length(id_Sym));% create a property for the id of the active slip system (parent)
grains.prop.activeSS=SS(changem(idSSactive,id_Sym,1:length(id_Sym)));% create a property for the active slip system (parent)

% 
% g_phase=grains(phaseName);
% colors=[ [0, 0.4470, 0.7410];...
%             [0.8500, 0.3250, 0.0980];...
%             [0.9290, 0.6940, 0.1250];...
%             [0.4940, 0.1840, 0.5560];...
%             [0.4660, 0.6740, 0.1880];...
%             [0.6350, 0.0780, 0.1840];...
%             [0.3010, 0.7450, 0.9330] ];
        
%     legend_str ={};proportionSS=[];
% 
%     for i=1:nSS
%         grains_active=g_phase(g_phase.prop.activeSSid_sym==i);
%         if ~isempty(grains_active)
%             grains_active.color=colors(i,:);
%             plot (grains_active);
%             hold on
%             proportionSS(i)=(length(grains_active)/length(g_phase))*100;
%             legend_str{i}= ['(',num2str(round(SS(i).n.hkl)) ,') [',num2str(round(SS(i).b.uvw)),'] (',num2str( round(proportionSS(i)) ) ,'%)'];
%         else
%             legend_str{i}={};
%         end  
%     end
%         legend_str = legend_str(~cellfun('isempty',legend_str)); % clear empty legend
%         legend(legend_str, 'fontSize',12)
%         hold off
%% Calculate mprime   
    gb_Gid=gb.grainId;
    SS_n=SS_sym(idSSactive).n;
    SS_b=SS_sym(idSSactive).b;

    n_in=SS_n(gb_Gid(:,1)); % slip plane normal of grains of one side of the boundary (miller idx)
    n_out=SS_n(gb_Gid(:,2)); % and the other side
    b_in=SS_b(gb_Gid(:,1)); %slip line of grains of one side of the boundary (miller idx)
    b_out=SS_b(gb_Gid(:,2)); % and the other side

    n_in=n_in./norm(n_in); % normalize vectors
    n_out=n_out./norm(n_out);
    b_in=b_in./norm(b_in);
    b_out=b_out./norm(b_out);

    n_in_rot=rotate(n_in, ori(gb_Gid(:,1))); % get orientation of slip plane normal in crystal frame
    n_out_rot=rotate(n_out,ori(gb_Gid(:,2)));
    b_in_rot=rotate(b_in,ori(gb_Gid(:,1))); % get orientation of slip direction in crystal frame
    b_out_rot=rotate(b_out,ori(gb_Gid(:,2)));

    mp=abs(dot (n_in_rot, n_out_rot) .* dot(b_in_rot, b_out_rot)); % m'
    
    rbv=norm(b_in_rot-b_out_rot); % magnitude of Residual Burger Vector