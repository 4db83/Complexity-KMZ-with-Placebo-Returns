function dir_name = set_dir(dir_name)
% call as: oos_dir = set_dir( './oos_eval/');

warning('off', 'MATLAB:MKDIR:DirectoryExists');
if ~exist(dir_name, 'dir')
  mkdir(dir_name);
end
warning('on', 'MATLAB:MKDIR:DirectoryExists');