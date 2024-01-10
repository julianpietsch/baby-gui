function exportDataToCSV( cExperiment, root_folder,name_root,detail_struct )
%exportDataToCSV( cExperiment, root_folder,name_root,detail_struct )
%exports data on the cellInf (the compiled cell info from the experiment) to a csv file for easier
%sharing.
%
%   cExperiment      -   cExperiment object
%   root_folder      -   (string) folder in which to save csv files
%   name_root        -   (string) first part of the final csv file names
%   detail_struct    -   (struct) structure encoding what to export and what the file should be called.
%                                 a structure array with the fields:
%                                       channel   -  (integer)channel to extract
%                                       field     -  (string)field of that channel to extract
%                                       name_tail -  (string)string to append to the name_root to make the whole .csv
%                                                    file name
%
% WARNING = if you have only extracted some of the channels you will have to put the channel number
% as the number of the cellInf structure entry you want to extract (so maybe just 1 if you only
% extracted 1 channel)
% precision is set to 16 dp (for the ordinate). Can set to higher if desired.


for filei = 1:length(detail_struct)
    detail_structi = detail_struct(filei);
    data_matrix = full(cExperiment.cellInf(detail_structi.channel).(detail_structi.field));
    dlmwrite(fullfile(root_folder,[name_root detail_structi.name_tail]),data_matrix,'delimiter',',','precision',16);
end

fprintf('finished writing data\n')

end

