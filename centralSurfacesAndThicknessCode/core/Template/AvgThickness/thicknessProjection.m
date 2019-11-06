function [EpProj]=thicknessProjection(sourceFile, targetFile, outputFileName)
% [EpProj]=thicknessProjection(sourceFile, targetFile, outputFileName)
% Project the thickness associated to a source mesh onto a target mesh
%
% Arguments:
%   sourceFile (str): path to the source mesh (with associated thickness
%       values)
%   targetFile (str): path to the target mesh
%   outputFileName (str): path where to stored the target mesh with
%       projected thickness values
%
% Returns:
%   EpProj (matrix): thickness values associated to each vertex of the
%       target mesh (after projection of the thickness from source to
%       target)

[Src.Vertices, Src.Faces, Src.Ep, Src.Colors, Src.TextureCoordinates] = VTKPolyDataReader(sourceFile);
[Tgt.Vertices, Tgt.Faces, Tgt.Ep, Tgt.Colors, Tgt.TextureCoordinates] = VTKPolyDataReader(targetFile);

Vert=cell(1,1);
Ep=cell(1,1);

Vert{1}=Src.Vertices'
Ep{1}=Src.Ep'

[EpProj]=projection_template_single(Vert,Ep,Tgt.Vertices,'wendland',4);

if strcmp(outputFileName,'None') ~= 1
  VTKPolyDataWriter(Tgt.Vertices, Tgt.Faces, EpProj, [Tgt.Ep Tgt.Ep Tgt.Ep], Tgt.Ep, outputFileName);
end
