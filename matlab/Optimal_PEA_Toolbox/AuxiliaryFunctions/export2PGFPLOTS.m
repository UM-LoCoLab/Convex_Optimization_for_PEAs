function export2PGFPLOTS( data2Export, dataLabels, nPoints, fileName)
%Generate .dat file from input data to generate plots in PGFPlots
%   export2PGFPLOTS( data2Export, DataLabels, nPoints, fileName) generates 
%   a .dat file in the current folder.
%
%       -- data2Export [type: struct]: each pair of fields in the structure
%       should include the (X,Y) data to be considered for PGFplots.
%       -- DataLabels [type: array of str]: each element of the array will
%       be used as labels to the data included in data2Export. The number 
%       of (X,Y) pairs in data2Export should match the number of labels.
%       -- nPoints [type: int]: number of data points to be included in the plot
%       -- fileName [type: str]: Name of the .dat file to save the results

fields = fieldnames(data2Export);

%-- Checking if data has been provided in pairs
if (rem(numel(fields), 2) == 0 && rem(numel(dataLabels), 2) == 0)
    numPairs = numel(fieldnames(data2Export))/2;
else
    error(['data2Export has to have an even number of elements to describe'...
        'the pairs of (X,Y) data points'])
end

%-- Checking if the number of DataLabels is consistent with number of datasets
if (numPairs) ~= numel(dataLabels)/2
    error('The amount of data labels does not match the data to export')
end

%-- Interpolating data and writting to .dat
for i = 1:numPairs
    ite = (i-1)*2+1;
    x = data2Export.(fields{ite});
    y = data2Export.(fields{ite+1});
    
    %-Enforcing monotonicity. Removing points from x that are not in order
    [~, indxX] = sort(x);
    [~, indxY] = sort(y);
    %-Substraction equal to 0 when the order aligns
    indxEqual = indxX - indxY;
    x = x(indxEqual == 0);
    y = y(indxEqual == 0);  
    %-Getting indexes for points that are different
    indxDifferentPoints = not(diff(x) == 0);
    x = x(indxDifferentPoints);
    y = y(indxDifferentPoints);
    %-Interpolating the x's and y's
    xInter = linspace(x(1), x(end), nPoints);
    yInter = interp1(x, y, xInter);
    
    %-- Step 2: Export .dat      
    %-- Exporting elongation and torque
    vec2Write = [xInter.', yInter.'];
    fid = fopen(sprintf('%s_%s.dat', fileName, dataLabels{ite+1} ),'w');
    fprintf(fid, sprintf('%s %s\n', dataLabels{ite}, dataLabels{ite+1}) );
    fprintf(fid,'%f %f\n', vec2Write.');
    fclose(fid);
end
