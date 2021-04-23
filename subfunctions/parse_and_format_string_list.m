function str2 = parse_and_format_string_list(str)
%input: a list of string separated by comas or smi colons
%output: a cell array of the substrings
%spaces are removed and ignored

%return empty output if empty input
if isempty(str)
    str2 = {[]};
    return
end

%remove spaces around commas
str = strrep(str,' ,', ',');
str = strrep(str,', ', ',');

%return empty output if empty input
if isempty(str)
    str2 = {[]};
    return
end

% remove leading/trailing spaces
while startsWith(str,' ')
    str = str(2:end);
end
while endsWith(str,' ')
    str = str(1:end-1);
end

% split along commas
str2 = strsplit(str,',');
end