function baseline_heightMatch = extract_nMLF_baseline(response)
    sti_repeat_number = 4;
    vol_rate = 1.2879;%FrameRate/(Zall+1);
    baseline_sec = 3;
    sti_dur = height(response)/sti_repeat_number;
    baseline_dur_idx = round(vol_rate*baseline_sec);
    basedex1 = (sti_dur-baseline_dur_idx):sti_dur;
    basedex2 = (sti_dur*2-baseline_dur_idx):sti_dur*2;
    basedex3 = (sti_dur*3-baseline_dur_idx):sti_dur*3;
    % basedex4 = (sti_dur*4-baseline_dur_idx):sti_dur*4;
    baseline = response(cat(2,basedex1,basedex2,basedex3),:);
    avgBaseline = median(baseline,1);
    baseline_heightMatch = {repelem(avgBaseline, height(response),1)};
end