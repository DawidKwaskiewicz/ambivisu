function [binauralAudio] = getBinaural(audio, fs, format)

ARIDataset = load('ReferenceHRTF.mat');
hrtfData = ARIDataset.hrtfData;
sourcePosition = ARIDataset.sourcePosition(:,[1,2]);

% Create a sphere with a distribution of points
nPoints = 24;   % number of points to pick
rng(0);         % seed randcom number generator
sphereAZ = 360*rand(1,nPoints);
sphereEL = rad2deg(acos(2*rand(1,nPoints)-1))-90;
pickedSphere = [sphereAZ' sphereEL'];

% Compare distributed points on the sphere to points from the HRTF dataset
pick = zeros(1, nPoints);
d = zeros(size(pickedSphere,1), size(sourcePosition,1));
for ii = 1:size(pickedSphere,1)
    for jj = 1:size(sourcePosition,1)
        % Calculate arc length
        d(ii,jj) = acos( ...
            sind(pickedSphere(ii,2))*sind(sourcePosition(jj,2)) + ...
            cosd(pickedSphere(ii,2))*cosd(sourcePosition(jj,2)) * ... 
            cosd(pickedSphere(ii,1) - sourcePosition(jj,1)));
    end
    [~,Idx] = sort(d(ii,:)); % Sort points
    pick(ii) = Idx(1);       % Pick the closest point
end

order = round(sqrt(size(audio, 2)) - 1);
devices = sourcePosition(pick,:)';
dmtrx = audioexample.ambisonics.ambidecodemtrx(order, devices);

FIR = cell(size(pickedSphere));
for ii = 1:length(pick)
    FIR{ii,1} = dsp.FrequencyDomainFIRFilter(hrtfData(:,pick(ii),1)');
    FIR{ii,2} = dsp.FrequencyDomainFIRFilter(hrtfData(:,pick(ii),2)');
end

desiredFs = 48e3;
audio = resample(audio,desiredFs,fs);

% if format == "fuma"
%     format2 = 'fuma-fuma';
% elseif format == "ambix"
%     format2 = 'acn-sn3d';
% else
%     format2 = 'error';
% end
audioFiltered = zeros(size(audio, 1),size(FIR,1),2);

audioDecoded = audioexample.ambisonics.ambidecode(audio, dmtrx, format);
for ii = 1:size(FIR,1)
    audioFiltered(:,ii,1) = step(FIR{ii,1}, audioDecoded(:,ii)); % Left
    audioFiltered(:,ii,2) = step(FIR{ii,2}, audioDecoded(:,ii)); % Right
end
audioOut = 10*squeeze(sum(audioFiltered,2));   % Sum at each ear 

binauralAudio = audioOut;

end