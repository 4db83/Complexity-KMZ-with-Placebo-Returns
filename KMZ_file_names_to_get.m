% --------------------------------------------------------------------------------------------------
% function to get existing file names to only compute those that don't already exist. Convenient
% when the cluster breaks down or runs out of time. 
% --------------------------------------------------------------------------------------------------
function list_of_seed_names = KMZ_file_names_to_get(OUTPUT_DIR, trnwin, stdize_Y, demean, nSim)

  % make function handle to easier create strings
  fsp = @(trnwin,stdize_Y,demean) ( strcat('maxP-12000-trnwin-', num2str(trnwin), ...
          '-gamma-2-stdize-', num2str(stdize_Y), '-demean-', num2str(demean), '-v2') );
  % string of file path
  save_path       = strcat( OUTPUT_DIR, fsp(trnwin,stdize_Y,demean) );
  sep('*'); 
  fprintf( save_path + ". "); 
  
  files_tmp       = dir([save_path]);
  files_listed    = vertcat(char(files_tmp.name));
  % GET FILE NAMES NUMBER ONLY AND CONVERT TO DOUBLE
  existing_seeds  = sort(str2double(strrep(strrep(cellstr(files_listed),'iSim',''),'.mat','')));
  seed_list_full  = (1:nSim)';
  file_list_diff  = setdiff(seed_list_full, existing_seeds)';
  % output
  list_of_seed_names = file_list_diff;

end 