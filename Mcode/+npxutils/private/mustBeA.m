function mustBeA(value, type)
%MUSTBEA Validate that an input is of a particular data type
%
% MUSTBEA(value, type)
%
% Validates that the input Value is of the specified Type or a
% subtype. If Value is not of Type, an error is raised. If Value is
% of Type, does nothing and returns.
%
% Value is the value to validates the type of. It may be anything. If
% you call it using a variable (as opposed to a longer expression),
% the variable name is included in any error messages.
%
% Type (char) is the name of the type that Value must be. A type
% name may be one of:
%   * A class, such as 'double', 'cell', or 'containers.Map'
%   * A Janklab pseudotype, such as 'cellstr' or 'numeric'
%
% Note: The cellstr pseudotype is nontrivial to check for, as it
% must call iscellstr() and check all cell contents. Avoid calling it in
% performance-critical code.
npxutils.internal.util.mustBeA(value, type);
end


