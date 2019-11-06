function save_thickness_figure(thickmap_path, thickfig_folder, thickfig_prefix, colourbar_min_str, colourbar_max_str)
% Save the thickness map visualisation from two different views as a
% .png image.
% This function is intended to be called by an external program
% (e.g. Python), hence the string conversions
%
% Arguments:
%    thickmap_path (str): path to the .mat structure containing the
%        central surface volume
%    thickfig_folder (str): path to the folder where to output .png
%        visualisation
%    thickfig_prefix (str): prefix for the name of the files that
%        will be created inside thickfig_folder
%    colourbar_min_str (str): min value for the colourbar. Should be as str
%        containing a float number
%    colourbar_max_str (str): max value for the colourbar. Should be as str
%        containing a float number


% Read central surface from file
thickmap = load(thickmap_path) ;

% convert string arguments to number
colourbar_min = str2num(colourbar_min_str) ;
colourbar_max = str2num(colourbar_max_str) ;

% save the thickness figure for thickness map
save_thicknessmap_figure(thickmap, thickfig_folder, thickfig_prefix, colourbar_min, colourbar_max) ;

end

