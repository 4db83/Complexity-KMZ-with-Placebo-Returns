function [] = start_parpool_with(NoWorkers)
% SET/CHECK HOW MANY CORES/LOGICAL PROCESSES: CTRL+SHIFT+ESC --> CPU % FEATURE('NUMCORES')
% Default = max No. Threads. 
% If NoWorkers is different to current, force stop.

c = parcluster;  
max_Numworkers = c.NumWorkers;

if nargin < 1 
  NoWorkers = max_Numworkers; 
end
% c.NumWorkers = NoWorkers; 
current_pool = gcp('nocreate');

% delete(gcp('nocreate')). % to close down
if isempty(current_pool) 
  c.NumWorkers = NoWorkers;
  p = c.parpool(NoWorkers); 
else
  % if (current_pool.NumWorkers < NoWorkers)||(force_new_pool == 1)
  if (current_pool.NumWorkers ~= NoWorkers)
    delete(current_pool);
    c.NumWorkers = NoWorkers; 
    p = c.parpool(NoWorkers);
  else
    fprintf('Parallel pool already has %d workers: \n', NoWorkers);
  end
end
