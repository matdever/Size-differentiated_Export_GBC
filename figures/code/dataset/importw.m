function wtotal = importw(filename)

%% Initialize variables.
delimiter = ',';
startRow = 1;
endRow = inf;

%% Format for each line of text:
formatSpec = '%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');


%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
wtotal = dataArray{:, 1};
end
