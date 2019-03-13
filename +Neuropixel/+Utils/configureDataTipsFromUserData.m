function configureDataTipsFromUserData(fig)
    if nargin == 0
        fig = gcf;
    end

    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@datatipfn);

    function txt = datatipfn(empt, event_obj)
        if ~isvalid(event_obj.Target), txt = ''; return; end
        ud = event_obj.Target.UserData;
        if ischar(ud)
            txt = ud;
        elseif isstruct(ud)
            x0 = getor(ud, 'xoffset', 0);
            xs = getor(ud, 'xscale', 1);
            xu = getor(ud, 'xunits', '');
            xn = getor(ud, 'xname', 'X');
            y0 = getor(ud, 'yoffset', 0);
            ys = getor(ud, 'yscale', 1);
            yu = getor(ud, 'yunits', '');
            yn = getor(ud, 'yname', 'Y');
            pos = event_obj.Position;
            xs = sprintf('%.3f %s', double(pos(1)-x0) * double(xs), xu);
            ys = sprintf('%.3f %s', double(pos(2)-y0) * double(ys), yu);
            ud.(xn) = xs;
            ud.(yn) = ys;
            
            flds = setdiff(fieldnames(ud), {'xoffset', 'yoffset', 'xscale', 'yscale', 'xunits', 'yunits', 'xname', 'yname'});
            txt = cell(numel(flds), 1);
            for iF = 1:numel(flds)
                fld = strrep(flds{iF}, '_', '\_');
                val = ud.(flds{iF});
                if isnumeric(val)
                    if isinteger(val)
                        val = sprintf('%d', val);
                    else
                        val = sprintf('%.3f', val);
                    end
                elseif islogical(val)
                    if val
                        val = 'true';
                    else
                        val = 'false';
                    end
                end
                txt{iF} = sprintf('\\color[rgb]{.15 .15 .15}\\rm%s \\color[rgb]{0 0.5 1}\\bf%s', fld, val);
            end
        else
            txt = '';
        end
        
        function v = getor(s, fld, def)
            if isfield(s, fld)
                v = s.(fld);
            else
                v = def;
            end
        end
    end
end