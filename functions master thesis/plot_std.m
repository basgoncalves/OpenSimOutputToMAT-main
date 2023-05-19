function plot_std(data, titles, alpha, color, smth, y_label)
[ax, pos,FirstCol,LastRow,LastCol]  = tight_subplotBG(size(data,2),0,[],[0.1 0.05],[0.1 0.05]);
F = linspace(1,100,size(data{1},2)); % x axis

for iDOF = 1:size(data,2)
    axes(ax(iDOF)); hold on; grid on;
    stdshade(data{iDOF}, [], color, F, smth)
    title(titles{iDOF}, 'Interpreter', 'none')
    if any(iDOF == FirstCol)
        ylabel(y_label)
    elseif any(iDOF == LastRow)
        xlabel('Gait Cycle [%]')
    end
end

tight_subplot_ticks (ax,LastRow,0)

mmfn
end