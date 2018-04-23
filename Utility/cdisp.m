function [] = cdisp(text,c)
    if(c), clc, end
    mydate = datestr(datetime);
    disp(['[' mydate ']']);
    if(iscell(text))
        cellfun(@disp,text);
    else
        disp(text);
    end
end