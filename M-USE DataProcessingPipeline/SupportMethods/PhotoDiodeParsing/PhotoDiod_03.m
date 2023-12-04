close all
% load('/Users/kia/Downloads/twh88_w_serial.mat')
sampleRate = 300;
Norm_win = 100;
Norm_winL = Norm_win*sampleRate;

NSamples = height(serialRecvData.Analog);%3E5; %Number of samples, here set to minimize the size
% if set to the lenght of the data gets the whole session data

photoL=serialRecvData.Analog.LightL(1:NSamples);
photoR=serialRecvData.Analog.LightR(1:NSamples);
timestamps = serialRecvData.Analog.SynchBoxTime(1:NSamples);
Event_ts = serialRecvData.EventCodes.SynchBoxTime;


EV_pts = find((Event_ts-Event_ts(1))/1E3<NSamples/sampleRate);
EVcodes = serialRecvData.EventCodes.CodeValue(EV_pts);
Event_ts=Event_ts(EV_pts);

% preprocessing the photodiode data
sampleL = length(photoL);
sampleR = length(photoR);
sizeL = size(photoL);
sizeR = size(photoR);

%ero-phase Butterworth bandpass filter with frequency cutoffs of 1Hz to 75Hz.
fL = bandpass(photoL,[1 75],sampleRate);
fR = bandpass(photoR,[1 75],sampleRate);


ConvWinL = ones(1,Norm_winL);
ConvWinR = gausswin(1,Norm_winL);
wrap_R = conv(abs(fR)',ConvWinR,'same')./conv(ones(1,sampleR),ConvWinR,'same');
wrap_L = conv(abs(fL)',ConvWinL,'same')./conv(ones(1,sampleL),ConvWinL,'same');

fL = fL./wrap_L';
fR = fR./wrap_R';

Norm_varL = (fL/max(abs(fL)));  
Norm_varR = (fR/max(abs(fR)));
Norm_varL = Norm_varL-nanmean(Norm_varL);
Norm_varR = Norm_varR-nanmean(Norm_varR);

logic_photoL = zeros(sizeL);
logic_photoR = zeros(sizeR);

[~,zero_Left] = findpeaks(-Norm_varL,'MinPeakDistance',6);
[~,one_Left] = findpeaks(Norm_varL,'MinPeakDistance',6);


logic_photoL = sign(Norm_varL);
logic_photoR = sign(Norm_varR);


zero_Right = find(diff(logic_photoR)<0)-2;
one_Right = find(diff(logic_photoR)>0)-2;

zero_Right(zero_Right>(sampleR-1))=[];
one_Right(one_Right>(sampleR-1))=[];


NumOffset = 2;
One_vec = zeros(1,sampleR+2*NumOffset);
Zero_vec = zeros(1,sampleR+2*NumOffset);
One_vec(one_Right+NumOffset)=1;
Zero_vec(zero_Right+NumOffset)=1;
inds = NumOffset+1:sampleR+NumOffset;
offset = [-NumOffset:NumOffset]';
inds_mat = repmat(inds,size(offset))+offset;
One_mat = One_vec(inds_mat);
Zero_mat = Zero_vec(inds_mat);

% while Go


% end



OxZ = sum(repmat(One_vec(NumOffset+1:end-NumOffset),size(offset)).* Zero_mat);
ZxO = sum(repmat(Zero_vec(NumOffset+1:end-NumOffset),size(offset)).* One_mat);

main_OneVec = One_vec(NumOffset+1:end-NumOffset);
main_ZeroVec = Zero_vec(NumOffset+1:end-NumOffset);


main_ZeroVec(find(ZxO==1)) = 0;
main_OneVec(find(OxZ==1)) = 0;
zero_Right = find(main_ZeroVec);
one_Right = find(main_OneVec);

OxZ = [];
ZxO = [];
main_ZeroVec = [];
main_OneVec = [];
One_mat = [];
Zero_mat = [];
inds_mat = [];


%% ploting the raw signals and detected transitions
startX = 20400;
xWidth = 800;

figure
subplot(2,2,1)
plot(Norm_varL,'LineWidth',1,'color',[.5 .5 .5])
hold on
plot(logic_photoL,'LineWidth',.5,'color','k')
%  plot(Rec_SignalL,'LineWidth',.5,'color','g')
set(gca,'ylim',[-2 2],'xlim',[startX startX + xWidth])
hold on
scatter(zero_Left,Norm_varL(zero_Left),20,'b')
scatter(one_Left,Norm_varL(one_Left),20,'r')
% scatter(adj_zeros,Norm_varL(adj_zeros),20,'filled','b')
% scatter(adj_ones,Norm_varL(adj_ones),20,'filled','r')


subplot(2,2,3)
hold on
plot(Norm_varR,'LineWidth',1,'color',[.5 .5 .5])
% plot(Rec_SignalR,'LineWidth',.5,'color','g') 
plot(logic_photoR,'LineWidth',.5,'color','k')
scatter(zero_Right,Norm_varR(zero_Right),20,'filled','r')
scatter(one_Right,Norm_varR(one_Right),20,'filled','b')
set(gca,'ylim',[-2 2],'xlim',[startX startX + xWidth])



subplot(2,2,2)
plot(photoL,'LineWidth',.5,'color','k')
hold on
scatter(zero_Left,photoL(zero_Left),20,'filled','b')
scatter(one_Left,photoL(one_Left),20,'filled','r')
set(gca,'xlim',[startX startX + xWidth])

subplot(2,2,4)
hold on
plot(photoR,'LineWidth',.5,'color','k')
scatter(zero_Right,photoR(zero_Right),20,'filled','r')
scatter(one_Right,photoR(one_Right),20,'filled','b')

set(gca,'xlim',[startX startX + xWidth])
set(gcf,'Position',[           1         458        1668         346])

%% Upsampling the left photodiode and making an ideal signal

sampleRate = 300;
fs2 = 1.5*1E3; % upsampling frequency

rateUp = fs2/sampleRate;
dt = 1/fs2;

Crr_offset = round(21*rateUp);
Crr_Win = fs2;

t_orig = (timestamps -timestamps(1));
t_r = t_orig(1):1/fs2:t_orig(end);
u_sample = length(t_r);
intp_Norm_varL = interp1(t_orig, Norm_varL, t_r');
intp_Norm_varR = interp1(t_orig, Norm_varR, t_r','pchip');
intp_photoR = interp1(t_orig, photoR, t_r','spline');

intp_logic_photoR = sign(intp_Norm_varR);

% finding upsampled transion points
[~,upsampled_zero_Left] = findpeaks(-intp_Norm_varL,'MinPeakDistance',6*rateUp);
[~,upsampled_one_Left] = findpeaks(intp_Norm_varL,'MinPeakDistance',6*rateUp);


% reconstructing the signal and detecting the transition points

t = [0:dt:(Crr_Win+Crr_offset-1)*dt]';
y_triangle = sin(2*pi*30*t);
signal_parts = 0:Crr_Win:u_sample;

if signal_parts(end)~= u_sample
   signal_parts = [signal_parts, u_sample];
end
Rec_SignalL =[];
for i = 1:length(signal_parts)-1
    p_Norm_vecL = intp_Norm_varL(signal_parts(i)+1:signal_parts(i+1));
    L_part = length(p_Norm_vecL);
    crr = zeros(1,Crr_offset);
    for j = 1:Crr_offset
       crr(j) = corr(y_triangle(j:L_part+j-1),p_Norm_vecL);
    end
    [~,ishift] =max(crr);
    Rec_SignalL = [Rec_SignalL;y_triangle(ishift:ishift+L_part-1)];
end

% adjusted upsample transition points according to the reconstructed ideal
% signal

[~,upsampled_adjusted_zeros] = findpeaks(-Rec_SignalL,'MinPeakDistance',8*rateUp,'MinPeakHeight',.4);
[~,upsampled_adjusted_ones] = findpeaks(Rec_SignalL,'MinPeakDistance',8*rateUp,'MinPeakHeight',.4);

%% Detect transition points on the right photodiode using the upsampled signal

upsamled_zero_Right = find(diff(intp_logic_photoR)<0)-round(2*rateUp);
upsampled_one_Right = find(diff(intp_logic_photoR)>0)-round(2*rateUp);


NumOffset = round(2*rateUp);
One_vec = zeros(1,u_sample+2*NumOffset);
Zero_vec = zeros(1,u_sample+2*NumOffset);
One_vec(upsampled_one_Right+NumOffset)=1;
Zero_vec(upsamled_zero_Right+NumOffset)=1;
inds = NumOffset+1:u_sample+NumOffset;
offset = [-NumOffset:NumOffset]';
inds_mat = repmat(inds,size(offset))+offset;
One_mat = One_vec(inds_mat);
Zero_mat = Zero_vec(inds_mat);

% while Go


% end



OxZ = sum(repmat(One_vec(NumOffset+1:end-NumOffset),size(offset)).* Zero_mat);
ZxO = sum(repmat(Zero_vec(NumOffset+1:end-NumOffset),size(offset)).* One_mat);

main_OneVec = One_vec(NumOffset+1:end-NumOffset);
main_ZeroVec = Zero_vec(NumOffset+1:end-NumOffset);



main_ZeroVec(find(ZxO)) = 0;
main_OneVec(find(OxZ)) = 0;
upsamled_zero_Right = find(main_ZeroVec);
upsampled_one_Right = find(main_OneVec);

OxZ = [];
ZxO = [];
main_ZeroVec = [];
main_OneVec = [];
One_mat = [];
Zero_mat = [];
inds_mat = [];


%% finding nearest adjusted upsample transition points to the non-adjusted transition points

% Find nearest matching sample
L_adj = sort([upsampled_adjusted_ones;upsampled_adjusted_zeros]);
L_raw = sort([upsampled_zero_Left;upsampled_one_Left] );
D = pdist2(L_adj, L_raw);

% For each point in A, find the nearest point in B
[min_distance, idx] = min(D, [], 2);

dt_adj = (t_r(L_adj)-t_r(L_raw(idx)));

%% reconstructing the right photodiode signal

% Rec_SignalR and Rec_SignalR2 are two reconstructed signals
% made by two different windows when cross correlating

% ishift and ishift2 shows the shift in number of sample for each
% cross-correlated window of the data, and iMax and iMax2 show maximum
% Xcorr. coeff. value of each segment of the data

RFv = [0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1];
NumCycle = 8 ;
orig_t = 0:1/60:NumCycle;
Rep_Sig = repmat(RFv',[round(length(orig_t)/length(RFv)),1])';
Crr_Win = 600;%round(length(Rep_Sig)/2);
Crr_offset = 600;

dt2 = 1/fs2; %delta T of upsampled signal
t_intp = [0:dt2:NumCycle]';
y_pChip = interp1(orig_t, [Rep_Sig,Rep_Sig(1)], t_intp,'nearest')-.5;
plot(y_pChip)

signal_parts = 0:Crr_Win:u_sample;

if signal_parts(end)~= u_sample
   signal_parts = [signal_parts, u_sample];
end
ishift =[];
iMax =[];
Rec_SignalR =[];
Rec_SignalR2 = [];

for i = 1:length(signal_parts)-1
    p_Norm_vecR = intp_logic_photoR(signal_parts(i)+1:signal_parts(i+1));
    L_part = length(p_Norm_vecR);
    crr = zeros(1,Crr_offset);
    for j = 1:Crr_offset
       crr(j) = corr(y_pChip(589+j:589+L_part+j-1),p_Norm_vecR,'type','Spearman');
    end
    [iMax(i),ishift(i)] =max(crr);
    
    Rec_SignalR = [Rec_SignalR;y_pChip(589+ishift(i):589+ishift(i)+L_part-1)];
    

end

Rec_SignalR2 = [];
signal_parts2 = 0:Crr_Win*5.5:u_sample;

if signal_parts2(end)~= u_sample
   signal_parts2 = [signal_parts2, u_sample];
end

for i = 1:length(signal_parts2)-1
    p_Norm_vecR = intp_logic_photoR(signal_parts2(i)+1:signal_parts2(i+1));
    L_part2 = length(p_Norm_vecR);
    crr = zeros(1,Crr_offset);
    for j = 1:Crr_offset
       crr(j) = corr(y_pChip(589+j:589+L_part2+j-1),p_Norm_vecR,'type','Spearman');
    end
    [iMax2(i),ishift2(i)] =max(crr);
    
    Rec_SignalR2 = [Rec_SignalR2;y_pChip(589+ishift2(i):589+ishift2(i)+L_part2-1)];
    

end


%%
% dt_adj = (t_r(L_adj)-t_r(L_raw(idx)));
d_shift = abs(diff(ishift'));
d_shiftR = min(d_shift,abs(d_shift-600));
disCR = (find(d_shiftR<1 ) )*Crr_Win;

%  disCR = (find( iMax2>.75 & iMax2<.85))*Crr_Win*5.5;

disC = disCR;
disC_t = t_r(disC(1:end-1));

iEx = 1;
jitterL = .5;

figure
subplot(2,2,1)
plot(t_r,intp_Norm_varL,'LineWidth',1,'color',[.5 .5 .5])
hold on
 plot(t_r,Rec_SignalL,'LineWidth',.5,'color','k')
% set(gca,'ylim',[-2 2],'xlim',[disC_t(iEx)-jitterL disC_t(iEx)+jitterL])
hold on
%  scatter(t_r(u_z_L),intp_Norm_varL(u_z_L),20,'b')
%  scatter(t_r(u_o_L),intp_Norm_varL(u_o_L),20,'r')
scatter(t_r(upsampled_adjusted_zeros),intp_Norm_varL(upsampled_adjusted_zeros),20,'filled','b')
scatter(t_r(upsampled_adjusted_ones),intp_Norm_varL(upsampled_adjusted_ones),10,'filled','r')

subplot(2,2,3)
hold on
% plot(t_r,intp_Norm_varR,'LineWidth',1,'color',[.5 .5 .5])
 plot(t_r,.5*Rec_SignalR+.5,'LineWidth',1,'color',[.2 .8 .2]) 
 plot(t_r,.5*Rec_SignalR2-.5,'LineWidth',1,'color','m') 
plot(t_r,intp_logic_photoR,'LineWidth',.5,'color','k')
plot(t_r,(intp_photoR-nanmean(intp_photoR))/100,'LineWidth',.5,'color','k')
scatter(t_r(upsamled_zero_Right),intp_Norm_varR(upsamled_zero_Right),20,'filled','r')
scatter(t_r(upsampled_one_Right),intp_Norm_varR(upsampled_one_Right),20,'filled','b')
set(gca,'ylim',[-2 2],'xlim',[disC_t(iEx)-jitterL disC_t(iEx)+jitterL])
set(gca,'ylim',[-2 2],'xlim',[disC_t(iEx)-jitterL disC_t(iEx)+jitterL])



%% find discrepencies
% 
% disC = L_adj(find((dt_adj)>=.012));
% disC_t = t_r(disC);
% 
% iEx = 55;
% jitterL = .5;
% figure
% subplot(2,2,1)
% plot(t_r,intp_Norm_varL,'LineWidth',1,'color',[.5 .5 .5])
% hold on
%  plot(t_r,Rec_SignalL,'LineWidth',.5,'color','k')
% set(gca,'ylim',[-2 2],'xlim',[disC_t(iEx)-jitterL disC_t(iEx)+jitterL])
% hold on
%  scatter(t_r(upsampled_zero_Left),intp_Norm_varL(upsampled_zero_Left),20,'b')
%  scatter(t_r(upsampled_one_Left),intp_Norm_varL(upsampled_one_Left),20,'r')
% scatter(t_r(upsampled_adjusted_zeros),intp_Norm_varL(upsampled_adjusted_zeros),20,'filled','b')
% scatter(t_r(upsampled_adjusted_ones),intp_Norm_varL(upsampled_adjusted_ones),10,'filled','r')
% 
% subplot(2,2,3)
% hold on
% % plot(t_r,intp_Norm_varR,'LineWidth',1,'color',[.5 .5 .5])
%  plot(t_r,Rec_SignalR,'LineWidth',.5,'color','g') 
% plot(t_r,intp_logic_photoR,'LineWidth',.5,'color','k')
% plot(t_r,(intp_photoR-nanmean(intp_photoR))/10,'LineWidth',.5,'color','k')
% scatter(t_r(upsamled_zero_Right),intp_Norm_varR(upsamled_zero_Right),20,'filled','r')
% scatter(t_r(upsampled_one_Right),intp_Norm_varR(upsampled_one_Right),20,'filled','b')
% set(gca,'ylim',[-2 2],'xlim',[disC_t(iEx)-jitterL disC_t(iEx)+jitterL])


%% plot histograms raw

% Find nearest matching sample

A = sort([zero_Right,one_Right])';
B = sort([zero_Left;one_Left] );
D = pdist2(A, B);

% For each point in A, find the nearest point in B
[min_distance, idx] = min(D, [], 2);

dt_adj = (timestamps(A)-timestamps(B(idx)))/1E3;

histogram(dt_adj,120,'FaceColor','b')

figure; 
subplot(3,1,1)
title dt-Frames-raw
histogram(diff(timestamps(sort([zero_Left;one_Left]))),80,'FaceColor','b')

% set(gca,'xlim',[0 100])

subplot(3,1,2)
title dt-Frames-adjusted
histogram(diff(timestamps(sort([one_Right,zero_Right]))),80,'FaceColor','r')
% set(gca,'xlim',[0 100])

subplot(3,1,3)
title dt-diff
histogram(dt_adj,120,'FaceColor','b')
% set(gca,'xlim',[0 100])

%% plot histograms

% Find nearest matching sample

A = sort([upsamled_zero_Right,upsampled_one_Right])';
B = sort([upsampled_adjusted_zeros;upsampled_adjusted_ones]);
D = pdist2(A, B);

% For each point in A, find the nearest point in B
[min_distance, idx] = min(D, [], 2);

dt_adj = t_r(A)-t_r(B(idx));



figure; 
subplot(4,1,1)
title dt-Frames-raw
histogram(diff(t_r(sort([upsampled_zero_Left;upsampled_one_Left]))),40,'FaceColor','b')
set(gca,'TickDir','out')
set(gca,'xlim',[0 .04])


subplot(4,1,2)
title dt-Frames-adjusted
histogram(diff(t_r(sort([upsampled_adjusted_ones;upsampled_adjusted_zeros]))),40,'FaceColor','b')
set(gca,'TickDir','out')
set(gca,'xlim',[0 0.04])

subplot(4,1,3)
title dt-Frames-adjusted
histogram(diff(t_r(sort([upsampled_one_Right,upsamled_zero_Right]))),80,'FaceColor','r')
set(gca,'TickDir','out')
set(gca,'xlim',[0 .15])

subplot(4,1,4)
title dt-diff
histogram(dt_adj,120,'FaceColor','b')
set(gca,'TickDir','out')
 set(gca,'xlim',[-.016 .016])


%%

figure
subplot(2,2,1)
plot(intp_Norm_varL,'LineWidth',1,'color',[.5 .5 .5])
hold on
 plot(Rec_SignalL,'LineWidth',.5,'color','k')
set(gca,'ylim',[-2 2],'xlim',[2E5*rateUp 2E5*rateUp+fs2])
hold on
 scatter(upsampled_zero_Left,intp_Norm_varL(upsampled_zero_Left),20,'b')
 scatter(upsampled_one_Left,intp_Norm_varL(upsampled_one_Left),20,'r')
scatter(upsampled_adjusted_zeros,intp_Norm_varL(upsampled_adjusted_zeros),20,'filled','b')
scatter(upsampled_adjusted_ones,intp_Norm_varL(upsampled_adjusted_ones),20,'filled','r')

subplot(2,2,3)
hold on
plot(intp_Norm_varR,'LineWidth',1,'color',[.5 .5 .5])
% plot(Rec_SignalR,'LineWidth',.5,'color','g') 
plot(intp_logic_photoR,'LineWidth',.5,'color','k')
scatter(upsamled_zero_Right,intp_Norm_varR(upsamled_zero_Right),20,'filled','r')
scatter(upsampled_one_Right,intp_Norm_varR(upsampled_one_Right),20,'filled','b')
set(gca,'ylim',[-2 2],'xlim',[2E5*rateUp 2E5*rateUp+fs2])



