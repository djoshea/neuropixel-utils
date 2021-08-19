function setLineOpacity(s, alpha)
% Set line opacity
%
% npxutils.internal.graphics.setLineOpacity(hLine, alpha)

for i = 1:numel(s)
    if ~verLessThan('matlab', '8.4')
        % use RGBA color specification
        if isvalid(s(i))
            s(i).Color(4) = alpha;
        end
    end
end

end