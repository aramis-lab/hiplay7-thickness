function computeAvgThickness_views(finalTemplateFileName, initialTemplateFileName, DataFileDir, InitSubjectsFileNames, TemplSubjectsFileNames, outputFilePrefix, testTitle)
% computeAvgThickness_views(finalTemplateFileName, initialTemplateFileName, DataFileDir, InitSubjectsFileNames, TemplSubjectsFileNames, outputFilePrefix, testTitle)
% Compute the average thickness (projected on a template) for controls
% and the average thickness (projected on the same template) for
% patients, as well as the difference in average thicknesses between the
% patients and the controls
%
% Arguments:
%   final_template_filename (str): path to the template obtained from
%       controls and a group of patients, using software Deformetrica
%   initial_template_filename (str): legacy, not used
%   DataFileDir (str): path where all the figures and data files output by
%       the current function will be stored
%   InitSubjectsFileNames (str): list of the surfaces for all the subjects
%       (controls and patients) used to produce final_template_filename. The
%       list is passed as a single string
%       "'[sub_1].vtk' '[sub_2].vtk' [...] '[sub_n].vtk'"
%   TemplSubjectsFileNames (str): list of the surfaces for all the
%       subjects after they have been deformed towards the template while
%       final_template_filename was computed with Deformetrica.
%       This is passed a a single string in a similar way to
%       InitSubjectsFilenames.
%   outputFilePrefix (str): prefix of the output data, mesh and image
%       files stored in DataFileDir
%   testTitle (str): legacy, not used
%
% Returns:
%   N/A

% Template data
[Template.Vertices, Template.Faces, Template.Ep, Colors, TextureCoordinates] = VTKPolyDataReader(finalTemplateFileName) ;

EpProj = [] ;
control_number = 0 ;
patient_number = 0 ;

idxPat = [] ;
idxControl = [] ;
idx = 1 ;


% Subjects data
InitSubjectsFileList=strsplit(InitSubjectsFileNames, ' ') ;
TemplSubjectsFileList=strsplit(TemplSubjectsFileNames, ' ') ;
nsubj=numel(TemplSubjectsFileList)

EpSuj=cell(1,nsubj) ;
for i=1:nsubj
   % Source of interpolation values
   TemplSubjectsFileList{i} ;
         
   TempEp = thicknessProjection(InitSubjectsFileList{i}, TemplSubjectsFileList{i}, strrep(TemplSubjectsFileList{i},'.vtk','_withThickness.vtk')) ;

   
   EpProj=[EpProj, TempEp];
   
   if numel(strfind(TemplSubjectsFileList{i}, 'patient')) > 0 
     idxPat=[idxPat, idx]
     patient_number=patient_number+1
   elseif numel(strfind(TemplSubjectsFileList{i}, 'control')) > 0
     idxControl=[idxControl, idx]
     control_number=control_number+1
   end
   
   idx=idx+1	
end

size(EpProj)
EpProj(1,:)

save(strcat(DataFileDir,'/',outputFilePrefix,'.mat'),'EpProj', 'idxPat', 'idxControl')

% FIGURE WITH ALL RESULTS

figh=figure
hold on

subplot(1,2,1, 'Position',[0 0.33 0.5 0.5]) %0 0.66 0.5 0.3
p=affiche_skel_3D_AF(Template.Faces, Template.Vertices, mean(EpProj(:,idxControl),2), 0, 2.5, [-90, 90]); % add template thickness to the mean!
view(25,15)
colormap jet ;
cmap=colormap;
axis off; lighting gouraud; camlight; 

%% Vue 2
h=subplot(1,2,2, 'Position',[0.5 0.33 0.5 0.5]) % 0.3 0.66 0.5 0.3
posh=get(h, 'pos')
posh(1)=posh(1) - 0.25
p=affiche_skel_3D_AF(Template.Faces, Template.Vertices, mean(EpProj(:,idxControl),2), 0, 2.5, [-90, 90]);
view(165,25)
axis off; lighting gouraud; camlight; 
colorbar

saveas(figh,strcat(DataFileDir,'/',outputFilePrefix,'_mean_controls_uniformscale_views.png'))
hold off

Col=mean(EpProj(:,idxControl),2);
VTKPolyDataWriter(Template.Vertices, Template.Faces, Col, [Col Col Col], Col, strcat(DataFileDir, '/', outputFilePrefix, '_mean_controls_uniformscale.vtk'));

figh=figure
colormap(cmap)
hold on

subplot(1,2,1, 'Position',[0 0.33 0.5 0.5]) %0 0.66 0.3 0.3 'Position',[0 0.33 0.5 0.3]
p=affiche_skel_3D_AF(Template.Faces, Template.Vertices, mean(EpProj(:,idxPat),2), 0, 2.5, [-90, 90]);
view(25,15)
axis off; lighting gouraud; camlight; 


h=subplot(1,2,2,'Position',[0.5 0.33 0.5 0.5]) %0.5 0.66 0.3 0.3 'Position',[0.3 0.33 0.5 0.3]
posh=get(h, 'pos')
posh(1)=posh(1) - 0.25
p=affiche_skel_3D_AF(Template.Faces,Template.Vertices, mean(EpProj(:,idxPat),2), 0, 2.5, [-90, 90]);
view(165,25)
axis off; lighting gouraud; camlight; 
colorbar

saveas(figh,strcat(DataFileDir,'/',outputFilePrefix,'_mean_TLE_uniformscale_views.png'))

hold off

Col=mean(EpProj(:,idxPat),2);
VTKPolyDataWriter(Template.Vertices, Template.Faces, Col, [Col Col Col], Col, strcat(DataFileDir, '/', outputFilePrefix, '_mean_TLE_uniformscale.vtk'));

figh=figure
colormap(cmap)
hold on
subplot(1,2,1, 'Position',[0 0.33 0.5 0.5]) % 0 0.33 0.3 0.3  'Position',[0 0 0.5 0.3]
p=affiche_skel_3D_AF(Template.Faces, Template.Vertices, mean(EpProj(:,idxControl),2)-mean(EpProj(:,idxPat),2), -0.75, 1.5, [-90, 90]);
view(25,15)
axis off; lighting gouraud; camlight; 


h=subplot(1,2,2, 'Position',[0.5 0.33 0.5 0.5]) % [0.5 0.33 0.3 0.3 'Position',[0.3 0 0.5 0.3]
posh=get(h, 'pos')
posh(1)=posh(1) - 0.5
p=affiche_skel_3D_AF(Template.Faces, Template.Vertices, mean(EpProj(:,idxControl),2)-mean(EpProj(:,idxPat),2), -0.75, 1.5, [-90, 90]);
view(165,25);
axis off; lighting gouraud; camlight; 
colorbar

saveas(figh,strcat(DataFileDir,'/',outputFilePrefix,'_mean_controls_patients_diff_uniformscale_views.png'))
hold off

Col=mean(EpProj(:,idxControl),2)-mean(EpProj(:,idxPat),2);
VTKPolyDataWriter(Template.Vertices, Template.Faces, Col, [Col Col Col], Col, strcat(DataFileDir,'/',outputFilePrefix, '_mean_diff_uniformscale.vtk'));

end

