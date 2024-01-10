function trap_prediction_image = generateTrapLocationPredictionImage(cTimelapse,timepoint,channel)
% function trap_prediction_image = generateTrapLocationPredictionImage(cTimelapse,timepoint,channel)
%
% generates trap prediction image by normalised cross correlation with
% trap template images in cTimelapse.trapTemplates.
%
% Takes 0.65*the maximum of the cross correlation as a threshold and picks all
% values below this in order, ruling out the area directly around each new
% trap Location.

image=cTimelapse.returnSingleTimepoint(timepoint,channel);
cTrap = cTimelapse.cTrapSize;
image=padarray(image,[cTrap.bb_height cTrap.bb_width],median(image(:)));

positiveExamples = cTimelapse.trapTemplates.positiveExamples;
posxcorr = cell(1,size(positiveExamples,3));
for i=1:size(positiveExamples,3)
  posxcorr{i} = normxcorr2(positiveExamples(:,:,i),image);
end
trap_prediction_image = mean(cat(3,posxcorr{:}),3);

if isfield(cTimelapse.trapTemplates,'negativeExamples')
    negativeExamples = cTimelapse.trapTemplates.negativeExamples;
    negxcorr = cell(1,size(negativeExamples,3));
    for i=1:size(negativeExamples,3)
      negxcorr{i} = normxcorr2(negativeExamples(:,:,i),image);
    end
    negxcorr = mean(cat(3,negxcorr{:}),3);
    trap_prediction_image = trap_prediction_image - negxcorr;
end

trap_prediction_image = ...
    trap_prediction_image(2*cTrap.bb_height+1:end-2*cTrap.bb_height,2*cTrap.bb_width+1:end-2*cTrap.bb_width);
f1 = fspecial('disk',1);
trap_prediction_image = imfilter((trap_prediction_image),f1,'same');
end
