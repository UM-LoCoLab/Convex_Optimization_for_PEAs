function  loadDataFromWeb()
%Load reference trajectories and datasets from a Dropbox folder.
%   The purpose of this file is to reduce the size of non-text files on the
%   git repository.
%   Edgar - August 2020

%Check if directory has been created
if ~exist('1_InputData_BiomechanicDatasets', 'dir')
    url = 'https://www.dropbox.com/sh/a3ter9ji7hqwt9p/AAA7-PcnWG4Wvc5yvhWhu6Sqa?dl=0';
    error(['No input data to run this code.\n'...
        'Please download the input data from the following Dropbox link'...
        ':\n\n%s'...
        '\n\nThen extract the files into the Optimal_SEA_Toolbox directory of this repo.'],...
        url);
end

end

