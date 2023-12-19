function discretizedFrames = GenerateDiscretizedFrameVector(FlipDetails)

discretizedFrames = [];

for iFlip = 1:height(FlipDetails)
    if ~ isnan(FlipDetails.NumFrames(iFlip))
        discretizedFrames = [discretizedFrames; ...
            [abs(zeros(FlipDetails.NumFrames(iFlip), 1) - FlipDetails.FlashPanelStatus(iFlip)) repmat(iFlip, FlipDetails.NumFrames(iFlip), 1)]];
    else
        discretizedFrames = [discretizedFrames; [NaN iFlip]];
    end
end





