function setMarkerOpacity(s, faceAlpha, edgeAlpha)
% stores information in UserData struct to cause saveFigure to render
% marker points as translucent when exporting to svg
    if nargin < 3
        edgeAlpha = faceAlpha;
    end

    for i = 1:length(s)
        % old version, simply tag it as translucent for saveFigure to pick
        % up during SVG authoring

        userdata = get(s(i),'UserData');
        userdata.svg.MarkerFaceAlpha = faceAlpha;
        userdata.svg.MarkerEdgeAlpha = edgeAlpha;
        set(s(i),'UserData', userdata);

        if ~verLessThan('matlab', '8.4')
            mh = s.MarkerHandle;
            if ~isa(mh, 'matlab.graphics.GraphicsPlaceholder')
                if ~isempty(mh.EdgeColorData)
                    mh.EdgeColorType = 'truecoloralpha';
                    mh.EdgeColorData(4) = uint8(edgeAlpha*255);
                end
                if ~isempty(mh.FaceColorData)
                    mh.FaceColorType = 'truecoloralpha';
                    mh.FaceColorData(4) = uint8(faceAlpha*255);
                end
            end

            % keep transparent
            addlistener(s(i),'MarkedClean',...
                @(ObjH, EventData) keepAlpha(ObjH, EventData, faceAlpha, edgeAlpha));
        end
    end
end

function keepAlpha(src, ~, faceAlpha, edgeAlpha)  
    mh = src.MarkerHandle;
    if ~isempty(mh.EdgeColorData)
        mh.EdgeColorType = 'truecoloralpha';
        mh.EdgeColorData(4) = uint8(edgeAlpha*255);
    end
    if ~isempty(mh.FaceColorData)
        mh.FaceColorType = 'truecoloralpha';
        mh.FaceColorData(4) = uint8(faceAlpha*255);
    end
end

