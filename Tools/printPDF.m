function printPDF(outputFn)

set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf, 'PaperPosition',[0 0 pos(3:4)],'PaperSize', pos(3:4));

if strcmpi(outputFn(end-3:end), '.pdf') == 0
    outputFn = strcat(outputFn, '.pdf');
end

print('-dpdf', outputFn);

