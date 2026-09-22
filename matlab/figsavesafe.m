function figsavesafe(fig, outFile)
% FIGSAVESAFE  Save a figure PNG robustly in headless batch mode.
% Tries saveas, then explicit painters print, then gives up gracefully
% (always closes the figure to avoid handle pileup).
    try
        saveas(fig, outFile);
    catch
        try
            print(fig, outFile, '-dpng', '-r150', '-painters');
        catch ME
            fprintf('figsavesafe: export failed for %s (%s)\n', outFile, ME.message);
        end
    end
    try, close(fig); catch, end
end
