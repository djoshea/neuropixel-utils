function out = todatetime(x, zonePolicy)
% Convert values to datetime array
%
% out = todatetime(x)
% out = todatetime(x, zonePolicy)
%
% TODATETIME converts an input array to a datetime array. It is an alternative
% to the datetime constructor, and is provided because MPT thinks Matlab got the
% datetime constructor behavior wrong: in the one-arg constructor, numeric inputs are
% interpreted as datevecs instead of datetimes; but we think that datenums are
% in more common use and should have been the default interpretation of numeric
% inputs. In addition, DATETIME does not support Java date/time types as inputs.
%
% TODATETIME also provides a zonePolicy argument for convenient setting
% of the TimeZone on the constructed datetimes.
%
% x may be any of:
%   - numeric, which are interpreted as datenums
%   - string, char or cellstr, which are parsed as datestrs
%   - datetime, which are left as is
%   - java.util.Date or java.util.Date[]
%   - java.time LocalDate, LocalDateTime, ZonedDateTime, OffsetDateTime, Instant, or
%       arrays of same.
%
% zonePolicy is a string indicating how the TimeZone on the constructed
% datetime array should be set. It may be one of:
%   - 'passthrough' (default) - Keep whatever zone was set on the input, if any
%   - 'unzoned' - no TimeZone is set
%   - 'utc' - TimeZone is set to UTC
%   - 'local' - TimeZone is set to the system's local time zone
%   - anything else - the name of a specific time zone to set on the output
%
% The way the zone policy is applied is:
%   * Inputs that do not have zone information are interpreted as being in the
%     zone specified by zonePolicy.
%   * Inputs that do have zone information are converted to the zone specified
%     by zonePolicy.
%
% If the input x is a datetime or zoned Java date and a zonePolicy is provided, 
% but the zonePolicy conflicts with the actual TimeZone set on x, an error is 
% raised.
%
% Returns a datetime array.
arguments
  x
  zonePolicy (1,1) string = 'passthrough'
end
out = npxutils.internal.util.todatetime(x, zonePolicy);
end