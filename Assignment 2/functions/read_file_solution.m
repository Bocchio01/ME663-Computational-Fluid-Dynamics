function [x, y] = read_file_solution(filepath)

fileID = fopen(filepath, 'r');
if fileID == -1
    error('Unable to open file');
end

x = [];
y = [];

try
    
    N_lines = sscanf(fgetl(fileID), 'ZONE I=%d, F=POINT');
    
    for i = 1:N_lines
        line = fgetl(fileID);
        parts = strsplit(line);
        x_val = str2double(parts{2});
        y_val = str2double(parts{3});
        x = [x; x_val];
        y = [y; y_val];
    end
    
    fclose(fileID);
    
catch
    fclose(fileID);
    error('Error reading file');
end

end
