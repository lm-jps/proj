function [data,emsg] = parse_json(string)
% [DATA,EMSG] = PARSE_JSON(string)
% This function parses a JSON string and returns a cell array with the
% parsed data. JSON objects are converted to structures and JSON arrays are
% converted to cell arrays.
% If EMSG is non-empty, an error occurred.  The string emsg documents
% the error.  On Matlabs that support "throw", an exception is thrown
% in this case instead.

% F. Glineur, 2009
% (inspired by the JSON parser by Joel Feenstra on MATLAB File Exchange
% (http://www.mathworks.com/matlabcentral/fileexchange/20565) but with 
% faster handling of strings)
% Turmon, 8/2010, add emsg support
% turmon, 11/2011, handle edge cases like zero-length responses

pos = 1;
len = length(string);
% String delimiters and escape characters are identified beforehand to improve speed
esc = regexp(string, '["\\]'); index_esc = 1; len_esc = length(esc);

emsg = ''; % signals OK
nc = next_char;
if isempty(nc),
  data = []; % needs to be valid
  error_pos('Outer level structure must be an object or an array, got nothing');
elseif nc == '{',
  data = parse_object;
elseif nc == '[',
  data = parse_array;
else,
  data = []; % needs to be valid
  error_pos('Outer level structure must be an object or an array');
end
return;

    function object = parse_object
        parse_char('{');
        object = [];
        if next_char ~= '}'
            while 1
                str = parse_string;
                if isempty(str)
                    error_pos(sprintf('Name of value at position %d cannot be empty', pos));
                end
                parse_char(':');
                val = parse_value;
                object.(valid_field(str)) = val;
                if next_char == '}'
                    break;
                end
                parse_char(',');
            end
        end
        parse_char('}');
    end

    function object = parse_array
        parse_char('[');
        object = cell(0, 1);
        if next_char ~= ']'
            while 1
                val = parse_value;
                object{end+1} = val;
                if next_char == ']'
                    break;
                end
                parse_char(',');
            end
        end
        parse_char(']');
    end

    function parse_char(c)
        skip_whitespace;
        if pos > len || string(pos) ~= c
            error_pos(sprintf('Expected %c at position %d', c, pos));
        else
            pos = pos + 1;
            skip_whitespace;
        end
    end

    function c = next_char
        skip_whitespace;
        if pos > len
            c = [];
        else
            c = string(pos);
        end        
    end
    
    function skip_whitespace
        while pos <= len && isspace(string(pos))
            pos = pos + 1;
        end
    end

     function str = parse_string
        % turmon fixed: check if we've overrun the string
        if pos > len,
            error_pos(sprintf('Error, expected string at position %d', pos));
        end
        if string(pos) ~= '"'
            error_pos(sprintf('String starting with " expected at position %d', pos));
        else
            pos = pos + 1;
        end
        str = '';
        while pos <= len
            while index_esc <= len_esc && esc(index_esc) < pos 
                index_esc = index_esc + 1;
            end
            if index_esc > len_esc
                str = [str string(pos:end)];
                pos = len + 1;
                break;
            else
                str = [str string(pos:esc(index_esc)-1)];
                pos = esc(index_esc);
            end
            switch string(pos)
                case '"' 
                    pos = pos + 1;
                    return;
                case '\'
                    if pos+1 > len
                        error_pos('End of JSON reached right after escape character');
                    end
                    pos = pos + 1;
                    switch string(pos)
                        case {'"' '\' '/'}
                            str(end+1) = string(pos);
                            pos = pos + 1;
                        case {'b' 'f' 'n' 'r' 't'}
                            str(end+1) = sprintf(['\' string(pos)]);
                            pos = pos + 1;
                        case 'u'
                            if pos+4 > len
                                error_pos('End of JSON reached in escaped unicode character');
                            end
                            str(end+1:end+6) = string(pos-1:pos+4);
                            pos = pos + 5;
                    end
                otherwise % should never happen
                    str(end+1) = string(pos);
                    pos = pos + 1;
            end
        end
        error_pos('End of JSON while expecting end of string');
    end

    function num = parse_number
        % turmon: changed arb. 20 below to 50 to be sure to catch entire number
        [num, one, err, delta] = sscanf(string(pos:min(len,pos+50)), '%f', 1); % TODO : compare with json(pos:end)
        if ~isempty(err)
            % turmon fixed bug here
            error_pos(sprintf('Error reading number at position %d', pos));
        end
        pos = pos + delta-1;
    end

    function val = parse_value
        % turmon fixed: check if we've overrun the string
        if pos > len,
            error_pos(sprintf('Error, expected value at position %d', pos));
        end
        switch(string(pos))
            case '"'
                val = parse_string;
                return;
            case '['
                val = parse_array;
                return;
            case '{'
                val = parse_object;
                return;
            case {'-','0','1','2','3','4','5','6','7','8','9'}
                val = parse_number;
                return;
            case 't'
                if pos+3 <= len && strcmpi(string(pos:pos+3), 'true')
                    val = true;
                    pos = pos + 4;
                    return;
                end
            case 'f'
                if pos+4 <= len && strcmpi(string(pos:pos+4), 'false')
                    val = false;
                    pos = pos + 5;
                    return;
                end
            case 'n'
                if pos+3 <= len && strcmpi(string(pos:pos+3), 'null')
                    val = [];
                    pos = pos + 4;
                    return;
                end
        end
        % turmon fixed bug here
        error_pos(sprintf('Value expected at position %d', pos));
    end

    % turmon -- altered this function for older matlab without throw
    % turmon -- altered to give more sane error message
    function error_pos(msg)
        if len < 2,
            snip = ['Whole JSON string = <', string, '>'];
        else,
          % extract the offending snippet of the input string
          % (this code needs the string to be at least 2 chars long)
          pad = 20; % padding before and after input snip
          dingbat = '<!err>';
          snip = strcat(string(max(pos-pad,1):min(pos,end)),...
                        dingbat,...
                        string(min(pos+1,end):min(pos+pad,end)));
          % indicate truncation
          if pos-pad > 1,
            snip = strcat('...', snip);
          end;
          if pos + pad > len,
            snip = strcat(snip, '...');
          end;
        end;
          
        % set up error message
        emsg = sprintf('%s.\nJSON context follows:\n%s', msg, snip);

        % For now, display it here to be sure it gets printed out once
        fprintf('parse_json: trouble.  Error: %s\n', emsg);
        
        % throw is a method of MException, and does not exist, in later matlabs
        if exist('throw') || exist('MException'),
          % keep old error-throwing behavior
          ME = MException('JSONparser:invalidFormat', emsg);
          throw(ME);
        else,
          % just return a nontrivial emsg, which the caller must check
          return;
        end;
    end

    function str = valid_field(str)   
    % From MATLAB doc: field names must begin with a letter, which may be
    % followed by any combination of letters, digits, and underscores.
    % Invalid characters will be converted to underscores, and the prefix
    % "alpha_" will be added if first character is not a letter.
        if ~isletter(str(1))
            str = ['alpha_' str];
        end
        str(~isletter(str) & ~('0' <= str & str <= '9')) = '_';   
    end

end
