function rgb = hsl2rgb(hsl)
% Convert Hue-Saturation-Luminance Color to Red-Green-Blue Color
%
% Usage:
%
%   RGB = npxutils.internal.graphics.hsl2rgb(HSL)
%
% Converts HSL, a M X 3 color matrix with values between 0 and 1
% into RGB, a M X 3 color matrix with values between 0 and 1
%
% See also: rgb2hsl, rgb2hsv, hsv2rgb

% Suresh E Joel, April 26,2003

narginchk(1, 1);

if max(max(hsl)) > 1 || min(min(hsl)) < 0
    error('HSL values have to be between 0 and 1');
end

for i = 1:size(hsl,1)
    if hsl(i,2) == 0 % when sat is 0
        % all values are same as luminance
        rgb(i,1:3) = hsl(i,3); %#ok<AGROW>
    end
    if hsl(i,3) < 0.5
        temp2 = hsl(i,3) * (1 + hsl(i,2));
    else
        temp2 = hsl(i,3) + hsl(i,2) - hsl(i,3) * hsl(i,2);
    end
    temp1 = 2*hsl(i,3) - temp2;
    temp3(1) = hsl(i,1) + 1/3;
    temp3(2) = hsl(i,1);
    temp3(3) = hsl(i,1) - 1/3;
    for j = 1:3
        if temp3(j) > 1
            temp3(j) = temp3(j) - 1;
        elseif temp3(j) < 0
            temp3(j) = temp3(j) + 1;
        end
        if 6 * temp3(j) < 1
            rgb(i,j) = temp1 + (temp2 - temp1) * 6 * temp3(j);
        elseif 2 * temp3(j) < 1
            rgb(i,j) = temp2;
        elseif 3 * temp3(j) < 2
            rgb(i,j) = temp1 + (temp2 - temp1) * (2/3 - temp3(j)) * 6;
        else
            rgb(i,j) = temp1;
        end
    end
end

% Sometimes the result is 1+eps instead of 1 or 0-eps instead of 0 ... so
% to get rid of this I am rounding to 5 decimal places)
rgb = round(rgb.*100000) ./ 100000; 
