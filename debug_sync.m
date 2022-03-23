good_cells = sp.cids(sp.cgs==2);
cellid = good_cells(29);
trigger_ms = patchCue_ts*1000;
spiket_ms = 1000*sp.st(sp.clu==cellid);
h=figure;
plot_timecourse('timestamp',spiket_ms,trigger_ms,trigger_ms-2000,trigger_ms+3000);
title(cellid);