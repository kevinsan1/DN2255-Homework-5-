function makeTable( firstColumn,secondColumn,thirdColumn,fourthColumn,fifthColumn,sixthColumn,tableName )
%makeTableLatex( firstColumn,secondColumn,thirdColumn,fourthColumn,tableName )
% Summary of this function goes here
%   Detailed explanation goes here
    currentD = cd;
	cd('/Users/kevin/SkyDrive/KTH Work/Period 3 2014/DN2255/Homework/5/Tables/')
%% If 6 columns
if nargin == 7
    FID = fopen(tableName, 'w');
    fprintf(FID, '\\begin{table}[h] \n');
    fprintf(FID, '   \\begin{center} \n');
    fprintf(FID, '      \\begin{tabular}{llllll}\\toprule \n');
    fprintf(FID, 'Variables\\\\ \n');
    fprintf(FID,'\\midrule \n');
    fprintf(FID,[...
        '$\\delta R(E)$ [$\\frac{\\text{mg}}{\\text{cm}^2}$]',...
        '& $\\delta E$ [keV]',...
        '& $\\delta R/R$ &',...
        '$\\delta R(E)$ [$\\frac{\\text{mg}}{\\text{cm}^2}$]',...
        '& $\\delta E$ [keV]',...
        '& $\\delta R/R$',...
        '\\\\ \\midrule \n']);
    for k=1:length(firstColumn)
        fprintf(FID, '%8.2f & %8.2f %8.2f & %8.2f %8.2f & %8.2f \\\\ ', firstColumn(k), secondColumn(k), thirdColumn(k), fourthColumn(k), fifthColumn(k),sixthColumn(k));
        if k==length(firstColumn)
            fprintf(FID, '\\bottomrule ');
        end
        fprintf(FID, '\n');
    end
    fprintf(FID, '      \\end{tabular} \n');
    fprintf(FID, '   \\end{center}\n');
    fprintf(FID,['\\caption{uncertainties in the range is calculated from',...
        'our extrapolation of the fit and',...
        ' Equation~\\eqref{eq:transformationEquations}',...
        ' $\\delta y = \\delta R/R$. The values obtained are',...
        ' around $1\\perc~$to$~2\\perc$, agreeing with Equation~\\eqref{eq:bergerUncertainty}.   } \n']);
    fprintf(FID,'\\label{tab:uncertaintyConversion} \n');
    fprintf(FID, '\\end{table} \n');
    fclose(FID);
else
    tableName = thirdColumn;
    t = num2cell([1,2,3,4]);
    [thirdColumn,fourthColumn,fifthColumn,sixthColumn] = deal(t{:});
    FID = fopen(tableName, 'w');
    fprintf(FID, '\\begin{table}[h] \n');
    fprintf(FID, '   \\begin{center} \n');
    fprintf(FID, '      \\begin{tabular}{llllll}\\toprule \n');
    fprintf(FID, '\\multicolumn{2}{c}{Variables}\\\\ \n');
    fprintf(FID,'\\midrule \n');
    for k=1:length(firstColumn)-1
        fprintf(FID, [firstColumn{k} ' & %4.4f \\\\ '], secondColumn(k));
        if k==length(firstColumn)-1
            fprintf(FID, [firstColumn{3} ' & %4.0f \\\\'],secondColumn(3));
            fprintf(FID, '\\bottomrule ');
        end
        fprintf(FID, '\n');
    end
    fprintf(FID, '      \\end{tabular} \n');
    fprintf(FID, '   \\end{center}\n');
    fprintf(FID,('\\caption{The values are used for Figure~\\ref{fig:Figures_correctedLevelSetMethod}.} \n'));
    fprintf(FID,'\\label{tab:myTable} \n');
    fprintf(FID, '\\end{table} \n');
    fclose(FID);
end
cd(currentD)
end

