function [ebsd,grains]=cleanGrains(ebsd,minCI,minANG,minGS,sF)
% Calculate grains, gb and plot IPF map
ebsd(ebsd.prop.ci<minCI).phaseId=-1;
ebsd(ebsd.prop.ci<minCI).phase=-1; 

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'threshold',minANG);

try
    notIndexed = grains('notIndexed');
    toRemove = notIndexed(notIndexed.grainSize ./ notIndexed.boundarySize<0.8);
    ebsd(toRemove) = [];
catch
end

ebsd(grains(grains.grainSize<minGS)) = [];

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'threshold',minANG);
grains=smooth(grains,sF);
end