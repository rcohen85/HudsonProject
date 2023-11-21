function [completeClicks, noise] = dHighres_click(p,hdr,bpDataHi)

% Tyack & Clark 2000 cite Au (1993) in Hearing by Whales & Dolphins, Au
% (ed.) stating that dolphins can distinguish clicks separated by as
% little as 205 us.

minGap_samples = ceil(p.mergeThr*hdr.fs/1e6);
energy = bpDataHi.^2;
candidatesRel = find(energy> (p.countThresh^2));
if length(candidatesRel)<1
    candidatesRel = [];
end
completeClicks = [];
noise = [];

if ~ isempty(candidatesRel)
    if p.saveNoise
        noise = dHR_get_noise(candidatesRel,length(energy),p,hdr);
    end
    [sStarts, sStops] = spDurations(candidatesRel, minGap_samples,length(energy));

    [c_starts,c_stops]= dHR_expand_region(p,hdr,...
            sStarts,sStops,energy,bpDataHi);
    
    completeClicks = [c_starts, c_stops];

end
