function writeGCFtoPDF(f, def_fn)
if nargin==1
    def_fn = 'Default.pdf';
end

[file, path] = uiputfile('*.pdf', 'Select File Name', def_fn);
if ~file
    return;
end

set(f, 'Renderer', 'Painters');
set(f, 'Units', 'inches');
figpos = get(f, 'Position');
set(f, 'Units', 'Normalized');
aspect = figpos(4)/figpos(3);

if (aspect<8/10.5)
    pos = [0 0 11 11*aspect];
else
    pos = [0 0 8.5/aspect 8.5];
end

set(f, 'PaperOrientation', 'Landscape');
set(f, 'PaperPosition', pos);
print(f, '-dpdf', fullfile(path,file));
