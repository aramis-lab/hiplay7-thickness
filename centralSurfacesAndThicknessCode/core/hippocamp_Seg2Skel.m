function hippocamp_Seg2Skel(dirData,dirResult,outputPrefix,sigma,export_VTK,export_matON,export_potON,export_thickON,tmpFile)
%hippocamp_Seg2Skel Central Surface and Thickness Map computation based on segmentations of the hippocampal strata pyramidale of the CA1-3 and the subiculum subregions.
%
%   Arguments:
%   - dirData: directory containing segmented images of hippocampal subregions in the BrainVisa IMA format. 
%        - the function expects each segmented subregion to be saved as an individual binary mask
%        - their respective filenames should contain a subregion identifier string: {'CA_GM','ca_gm','CA_SP','ca_sp','subiculum','SUBICULUM','CA_WM','ca_wm','sub_WM','SUB_WM','DG','dg','alveus','ALVEUS','fimbria','FIMBRIA'} . 
%   - dirResult: output directory 
%   - outputPrefix: output file prefix
%   - sigma(optional): size of the kernel (influences the regularity of the vector field). Default value=10
%   - export_VTK(optional): value=1 (default) to save the central surfaces in VTK polydata format; value=0, otherwise
%   - bool export_matOn(optional):  value=1 (default) to save the environment variables in a .mat file; value=0, otherwise
%   - bool export_potOn(optional):  value=1 (default) to save the potential as a brainvisa .ima file; value=0, otherwise
%   - bool export_thickOn(optional):  value=1 (default) to save the thickness maps as a brainvisa .ima file; value=0, otherwise
%   - tmpFile(optional) : value=0 (default), or =1 to save a temporary file containing the minimum and maximum thickness values


% initialise data to default values if not provided by user
if ~exist('sigma')
     sigma = 10;
end  
if ~exist('tmpFile')
   tmpFile='';
end  
if ~exist('export_thickON')
   export_thickON = 0;
end  
if ~exist('export_potON')
   export_potON = 0;
end  
if ~exist('export_matON')
   export_matON = 0;
end  
if ~exist('export_VTK')
    export_VTK = 0;
end  


% Read data from input folder
%-- Vol: hippocampal ribbon of gray matter
%-- CA_GM: gray matter of CA
%-- subiculum: gray matter of subiculum
%-- Eps_init: defines the orientation of the normal vector to the ribbon boundaries
%-- metadata: voxel size, volume dimensions
[Vol, CA_GM, subiculum, Eps_init, metaData] = generateInitialEpsilonFromRibbonBoundaries(dirData, dirData, 1);

opt.sigma = sigma ;
min_vox_size=min(metaData.vox_size) ;
for axis_index = [1:3]
    opt.aniso(axis_index) = min_vox_size/metaData.vox_size(axis_index) ;
end

% Compute 3D vector field within Vol
%-- V: vector field
%-- T1: length at each x of Vol of the streamline from inner boundary to x
%-- T2: length at each x of Vol of the streamline from x to outer boundary
%-- B1: all streamlines originating at any point within the domain and following
%       the vector field
%-- B2: all streamlines originating at any point within the domain and following
%       the opposite vector field
[V,T1,T2,B1,B2] = calcul_champ_3D(Vol.mat,Eps_init.mat,CA_GM.mat,subiculum.mat,opt);

%% Compute central surface
[Skel,EpVox,SkelFace,SkelVert]=calcule_skel_3D(Vol.mat,T1,T2,opt.aniso);

mm_origin = (metaData.origin-1.0).*metaData.vox_size ;
%SkelVert = repmat(mm_origin, [size(SkelVert,1),1])+SkelVert(:, [2,1,3]) * opt.aniso*metaData.vox_size ;
SkelVert = repmat(mm_origin, [size(SkelVert,1),1])+SkelVert(:, [2,1,3])*min_vox_size ;

% adjust thickness according to voxel size
SkelThick = EpVox.*min_vox_size;

%% Export variables as .mat file
if export_matON
    file_mat = strcat(dirResult,'/',outputPrefix);
    option.sigma = sigma;
    option.vox_size = metaData.vox_size;
    %dist_map.out = T1;
    %dist_map.in = T2;
    Vect_field = V;
    fprintf('save matlab data as %s \n',strcat(outputPrefix,'.mat'));
    save(file_mat, 'Skel','SkelVert','SkelFace','SkelThick','Vect_field','option','Vol','Eps_init','T1','T2','B1','B2','metaData');
end

%% Export surface as brainvisa .mesh file
if dirResult(end) ~= '/'
    dirResult = strcat(dirResult,'/');
end

file_mesh = strcat(dirResult,outputPrefix,'.mesh');
fprintf('export folder: %s \n', dirResult);
fprintf('export surface as %s.%s...', outputPrefix,'mesh');
exportMesh(file_mesh,SkelVert,SkelFace,SkelThick);
fprintf(' DONE \n');

%% Optional exports
% Export surface as .VTK
if export_VTK==1
    file_mesh = strcat(dirResult,outputPrefix,'.vtk');
    fprintf('export folder: %s \n', dirResult);
    fprintf('export surface as %s.%s...', outputPrefix,'vtk');
    VTKPolyDataWriter(SkelVert, SkelFace, SkelThick,[SkelThick SkelThick SkelThick], SkelThick, file_mesh); % VTKPolyDataWriter(SkelVert, SkelFace, EpVox,[EpVox EpVox EpVox],EpVox, file_mesh); 
end

% Export thickness potentials and maps
[m,n,l] = size(Vol.mat)
size(T1)
metaData.origin
% Export thickness potentials in the range of (-1 1)
if export_potON   
    tempIma.mat=zeros(metaData.dim_mat);
    tempIma.mat(metaData.origin(1):metaData.origin(1)+m-1,metaData.origin(2):metaData.origin(2)+n-1,metaData.origin(3):metaData.origin(3)+l-1) = ((T1-T2)>0) - ((T1-T2)<0);
    tempIma.dim_mat=size(tempIma.mat);
    tempIma.vox_size=metaData.vox_size;
    name_file=strcat(dirResult,outputPrefix,'_pot.ima')
    fprintf('export potential as %s \n', name_file);
    save_ima(name_file,tempIma,'DOUBLE','b')
end

% Export thickness map as brainvisa .ima file
if export_thickON    
    tempIma.mat=zeros(metaData.dim_mat);   
    tempIma.mat(metaData.origin(1):metaData.origin(1)+m-1,metaData.origin(2):metaData.origin(2)+n-1,metaData.origin(3):metaData.origin(3)+l-1) = T1 + T2;
    tempIma.dim_mat=size(tempIma.mat);
    tempIma.vox_size=metaData.vox_size;
    name_file=strcat(dirResult,outputPrefix,'_thick.ima')
    fprintf('export thickness map as %s \n', name_file);
    save_ima(name_file,tempIma,'DOUBLE','b')
end

% Export skeleton as brainvisa .ima file
if export_thickON    
    %-- compute skeleton in cropped volume
    skelmat = zeros(size(Vol.mat)) ;
    vert_number = size(Skel,1) ;
    for vert_index=[1:vert_number]
        skel_xyz = Skel(vert_index,:) ;
        skel_x = round(skel_xyz(1)) ;
        skel_y = round(skel_xyz(2)) ;
        skel_z = round(skel_xyz(3)) ;
        skelmat(skel_x, skel_y, skel_z) = 1 ;
    end
    %-- un-crop
    tempIma.mat=zeros(metaData.dim_mat);   
    tempIma.mat(metaData.origin(1):metaData.origin(1)+m-1,metaData.origin(2):metaData.origin(2)+n-1,metaData.origin(3):metaData.origin(3)+l-1) = skelmat;
    tempIma.dim_mat=size(tempIma.mat);
    tempIma.vox_size=metaData.vox_size;
    name_file=strcat(dirResult,outputPrefix,'_skel.ima')
    fprintf('export skeleton volume as %s \n', name_file);
    save_ima(name_file,tempIma,'DOUBLE','b')
end

% Export interface data
tmpFile=strcat(dirResult,tmpFile)
if numel(tmpFile) > 0
fid = fopen(tmpFile,'w');              	    		  % Open the file
if fid ~= -1
  fprintf(fid,'%f|%f\n', min(SkelThick), max(SkelThick));        % Print
  fclose(fid);                      	    		  % Close the file
end
end

end

