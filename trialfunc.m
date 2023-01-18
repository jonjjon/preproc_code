function [trl, event] = trialfunc(cfg);
    trl = [];
    
    % read the header information and the events from the data
    hdr   = ft_read_header(cfg.dataset);
    event = ft_read_event(cfg.dataset);

    % compare event value
    onevent = zeros(1,length(event));
    for i1 = 1:length(cfg.trialdef.eventvalue)
        onevent(find(strcmp(cfg.trialdef.eventvalue{i1},{event.value}))) = i1;
    end
    
    % determine the number of samples before and after the trigger
    ievent = find(onevent);
    for i1 = 1:length(ievent),
        begsample = event(ievent(i1)).sample + cfg.trialdef.prestim  * hdr.Fs;
        endsample = event(ievent(i1)).sample + cfg.trialdef.poststim * hdr.Fs-1; % adjustment for including 0
        offset    = cfg.trialdef.prestim  * hdr.Fs;
        trigger   = onevent(ievent(i1));
        trl(end+1,:) = [round([begsample endsample offset])  trigger];
    end
    
    return;
