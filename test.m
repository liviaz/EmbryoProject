% test read file for javs

f = fopen('C:\Users\Livia\Desktop\Sample Read.txt', 'r');
proceed = true;
numbers = [];

% read until we get to "Values are ..." line
while(~feof(f) && proceed)
    
    tline = fgetl(f);
    
    if length(tline) > 10 && isequal(tline(1:10), 'Values are')
        proceed = false;
    end
end

% read one more blank line, then start bring in numbers
tline = fgetl(f);

% start reading numbers
while (~feof(f))
    
    % parse string and add to numbers array
    [a, count] = fscanf(f, '%s', 1);
    eLoc = strfind(a, 'E');
    
    if eLoc > 2
        aNum = str2double(a(1:(eLoc - 1)));
        aExp = str2double(a((eLoc + 1):end));
        numbers = [numbers (aNum * 10^(aExp))];
    end
end

fclose(f);