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
if lower(zonePolicy) == "UTC"
  zonePolicy = "utc";
end

% TODO: Add support for Joda Time for legacy code?

% Convert values

if isa(x, 'datetime')
  out = x;
elseif isnumeric(x)
  out = datetime(x, 'ConvertFrom','datenum');
elseif isstring(x) || ischar(x) || iscellstr(x)
  % TODO: Replace this with our own try/catch-less parsing logic
  out = datetime(x);
elseif isjava(x)
  out = convertJavaDateToDatetime(x, zonePolicy);
else
  % Anything else, let method overrides or datetime's own behavior sort it out
  out = datetime(x);
end

% Apply TimeZone policy

if zonePolicy == "passthrough"
  % NOP; leave whatever zone was set on input
elseif zonePolicy == "unzoned"
  out.TimeZone = [];
elseif zonePolicy == "utc"
  out.TimeZone = 'UTC';
elseif zonePolicy == "local"
  out.TimeZone = datetime.SystemTimeZone;
else
  % Looks like a specific TimeZone was requested
  out.TimeZone = zonePolicy;
end
  
end

function out = convertJavaDateToDatetime(j, zonePolicy)
% TODO: Once this code has been validated, optimize it by converting it to
% use toInstant() and toEpochSecond() on java.time datetime values.
nanoScale = 10^9;
if isa(j, 'java.util.Date')
  out = datetime(j.getTime / 1000, 'ConvertFrom','posixtime');
  tzoff = j.getTimezoneOffset;
  if tzoff ~= 0
    if tzoff > 0
      sign = '-';
    else
      sign = '+';
    end
    absoff = abs(tzoff);
    offhours = floor(absoff / 60);
    offminutes = rem(absoff, 60);
    tzoffStr = sprintf('%s%02d:%02d', sign, offhours, offminutes);
    out.TimeZone = 'UTC';
    if zonePolicy == "utc"
      % Leave as is (optimization)
    else
      % Make the zone on the input date visible
      out.TimeZone = tzoffStr;
    end
  end  
elseif isa(j, 'java.time.Instant')
  sec = j.getEpochSecond;
  plusNanos = j.getNano;
  posixTime = double(sec) + (double(plusNanos) / nanoScale);
  out = datetime(posixTime, 'ConvertFrom','posixtime');
  % Instants are in UTC by definition
  out.TimeZone = 'UTC';
elseif isa(j, 'java.time.LocalDate')
  out = convertJavaDateToDatetime(j.atStartOfDay);
elseif isa(j, 'java.time.LocalDateTime')
  secondsWithNanos = double(j.getSecond) + (double(j.getNano) / nanoScale);
  out = datetime(j.getYear, j.getMonthValue, j.getDayOfMonth, j.getHour, j.getMinute, secondsWithNanos);
elseif isa(j, 'java.time.ZonedDateTime')
  secondsWithNanos = double(j.getSecond) + (double(j.getNano) / nanoScale);
  out = datetime(j.getYear, j.getMonthValue, j.getDayOfMonth, j.getHour, j.getMinute, secondsWithNanos);
  out.TimeZone = char(j.getZone.getId);
elseif isa(j, 'java.time.OffsetDateTime')
  secondsWithNanos = double(j.getSecond) + (double(j.getNano) / nanoScale);
  out = datetime(j.getYear, j.getMonthValue, j.getDayOfMonth, j.getHour, j.getMinute, secondsWithNanos);
  out.TimeZone = char(j.getOffset.getId);  
elseif isa(j, 'java.util.Date[]') || isa(j, 'java.time.Instant[]') || isa(j, 'java.time.LocalDate[]') ...
    || isa(j, 'java.time.LocalDateTime[]') || isa(j, 'java.time.ZonedDateTime[]') ...
    || isa(j, 'java.time.OffsetDateTime[]')
  % TODO: Push this loop down in to Java for efficiency?
  out = repmat(NaT, [1 numel(j)]);
  zone = [];
  for i = 1:numel(j)
    jEl = j(i);
    if isempty(jEl)
      % It's a Java null; leave as NaT
    else
      dt = convertJavaDateToDatetime(jEl);
      if isempty(zone)
        zone = dt.TimeZone;
        out.TimeZone = dt.TimeZone;
      else
        if ~isequal(dt.TimeZone, zone)
          error(['Input java.util.Date[] array has inconsistent time zones on its ' ...
            'elements; cannot convert to datetime']);
        end
      end
      out(i) = dt;
    end
  end
elseif isa(j, 'java.util.Collection')
  out = repmat(NaT, [1 j.size]);
  zone = [];
  it = j.iterator;
  i = 1;
  while it.hasNext
    jEl = it.next;
    if isempty(jEl)
      % It's a Java null; leave as NaT
    else
      dt = convertJavaDateToDatetime(jEl);
      if isempty(zone)
        zone = dt.TimeZone;
        out.TimeZone = dt.TimeZone;
      else
        if ~isequal(dt.TimeZone, zone)
          error(['Input Collection has inconsistent time zones on its ' ...
            'elements; cannot convert to datetime']);
        end
      end
      out(i) = dt;
    end
    i = i + 1;
  end
else
  error('Conversion to datetime from Java type %s is not supported', class(j));
end

end