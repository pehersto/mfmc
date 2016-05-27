%
% See COPYRIGHT for license information 
%
function MYDEBUG( msg, cond )
%MYDEBUG Prints debug message if cond satisfied
% In
%   msg     ...     message
%   cond    ...     msg printed if cond evaluates to true/1

if(cond)
    disp(msg);
end


end

