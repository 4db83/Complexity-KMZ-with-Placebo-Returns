folderName = 'KMZ-files_stdize_RFFs_0';           % folder you want to zip
zipName    = 'KMZ-files_stdize_RFFs_0.zip';       % output zip file

% delete any old archive first
if exist(zipName, 'file')
    delete(zipName);
end

% create fresh zip with "store only" (-0)
cmd = sprintf('zip -0 -r %s %s', zipName, folderName);
tic; status = system(cmd); toc;