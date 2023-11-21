function [startsSec,stopsSec] = dST_choose_segments_raw(hdr)

dnum2sec = 60*60*24;
starts = hdr.raw.dnumStart;
stops = hdr.raw.dnumEnd;
startsSec = (starts - starts(1)).*dnum2sec;
stopsSec = (stops - starts(1)).*dnum2sec;

