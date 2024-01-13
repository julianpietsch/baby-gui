function trackTrapsThroughTime(cTimelapse,timepoints,preserve_trap_info,image_reset_time,huge_move,image_frac)
% trackTrapsThroughTime(cTimelapse,timepoints,preserve_trap_info)
% 
% cTimelapse            -  an object of the babyTimelapse class.
% timepoints            -  (optional) an array of timepoint to track. Default is
%                          cTimelapse.timepointsToProcess. Does not have to be
%                          consecutive.
% preserve_trap_info    -  (optional) boolean: false to replace the
%                          trapInfo and trapLocations structures of
%                          cTimelapse.cTimepoint(timepoints) with the blank
%                          trapInfo structure. If true, only fills in empty
%                          timepoints.
%                          Default false - so replaces by default.
% image_reset_time      -  time after which image used for cross
%                          correlation is reset. default 80.
%
% 
% uses cross correlation to track the location of the traps at each
% timepoint. Start with timpoints(1) and cross correlates with the central
% region of this image until either the traps have moved half a trap width
% or image_reset_time timpoints have elapsed. Then resets to the current
% timepoint and uses cross correlation of images with this timpoint to
% estimate image drift. This image drift is added to the location of the
% traps at timepoints(1) to get the location of the traps at each timpoint.
%
% The locations of the traps at each timepoint in each image are stored in
% the trapLocation field of the cTimepoint structure. They are stored as a
% structure array with fields xcenter and ycenter: one element of the
% structure array for each trap.
% preexisting trapLocations are overwritten unless preserve_trap_info is
% true. trapInfo is similarly instantiated/replaced at all timepoints
% unless preserve_trap_info is true, in which case only those it is
% instantiated only when it is empty or absent.
%
% Uses cTimelapse.trapTemplateChannel channel for cross correlation.
%
% If traps are not present (cTimelapse.trapsPresent = false) then no drift
% is measured and it just instantiates the trapInfo structure, setting the
% cTimpoint structure to have just one trapInfo structure and a
% trapLocation of [0,0]. preserve_trap_info the trapInfo in
% just the same way as if traps are present.
%
% WARNING: if the drift is less than 1 pixel it will not be detected, so if
% it is drifting very slowly the value of image_reset_time may have to be
% increased.
%

if nargin<2 || isempty(timepoints)
    timepoints=cTimelapse.timepointsToProcess;
end

if nargin<3 || isempty(preserve_trap_info)
    preserve_trap_info=false;
end

if nargin<4 || isempty(image_reset_time)
    image_reset_time = 80;
end

if nargin<5 || isempty(huge_move)
    if all(isfield(cTimelapse.cTrapSize,{'bb_width','bb_height'}))
        huge_move = [cTimelapse.cTrapSize.bb_width/2,...
            cTimelapse.cTrapSize.bb_height/2];
    else
        huge_move = [cTimelapse.imSize(2)/4,cTimelapse.imSize(1)/4];
    end
end

if nargin<6 || isempty(image_frac)
    image_frac = 0.8;
end

trapInfo_struct = cTimelapse.createTrapInfoTemplate;

if cTimelapse.trapsPresent
    
    timepoint = timepoints(1);
    % take central region of first timepoint for registration.
    reg_im=cTimelapse.returnSingleTimepoint(timepoint,cTimelapse.trapTemplateChannel);
    bb=floor(size(reg_im)*(1-image_frac));
    accum_col=0;
    accum_row=0;
    
    reg_im=double(reg_im);
    reg_im=reg_im(bb(1):end-bb(1),bb(2):end-bb(2));
    reg_im_fft=fft2(reg_im/mean(reg_im(:)));
    timepoint_reg=timepoint;
        
    %make trapInfo_struct the size of the the number of traps at size
    %cTimelapse.cTimepoint(timpoint(1))
    trapInfo_struct(1:length(cTimelapse.cTimepoint(timepoint).trapLocations)) = trapInfo_struct;
    
    if ~preserve_trap_info || ~isfield(cTimelapse.cTimepoint(timepoint),'trapInfo') || isempty(cTimelapse.cTimepoint(timepoint).trapInfo)
        cTimelapse.cTimepoint(timepoint).trapInfo = trapInfo_struct;
    end
    
    for i=2:length(timepoints)
        timepoint=timepoints(i);
        
        % Trigger the TimepointChanged event for babyLogging
        babyLogging.changeTimepoint(cTimelapse,timepoint);
        
        new_im=cTimelapse.returnSingleTimepoint(timepoint,cTimelapse.trapTemplateChannel);
        
        new_im=new_im(bb(1):end-bb(1),bb(2):end-bb(2));
        new_im_fft2 = fft2(new_im/mean(new_im(:)));
        
        [row_dif, col_dif] = dftregistration(reg_im_fft,new_im_fft2);
                
        %If a huge move is deteected at a single timepoint it is taken to
        %be inaccurate and the correction from the previous timepoint used
        %(this might be common if there is a focus loss for example).
        if abs(col_dif-accum_col)>huge_move(1)
            col_dif=accum_col;
        end
        if abs(row_dif-accum_row)>huge_move(2)
            row_dif=accum_row;
        end
        
        accum_col = col_dif;
        accum_row = row_dif;
        
        
        xloc=[cTimelapse.cTimepoint(timepoint_reg).trapLocations(:).xcenter]-col_dif;
        yloc=[cTimelapse.cTimepoint(timepoint_reg).trapLocations(:).ycenter]-row_dif;
        
        %keep traps located on the image.
        xloc(xloc<1) = 1;
        xloc(xloc>cTimelapse.imSize(2)) = cTimelapse.imSize(2);
        yloc(yloc<1) = 1;
        yloc(yloc>cTimelapse.imSize(1)) = cTimelapse.imSize(1);
        
        % ensure pre existing trap locations will not be
        % overwritten by new trapLocations if preserve_trap_info is true.
        % Important for continuous segmentation.
        if ~preserve_trap_info || ~isfield(cTimelapse.cTimepoint(timepoint),'trapLocations') || isempty(cTimelapse.cTimepoint(timepoint).trapLocations)
            cTimelapse.cTimepoint(timepoint).trapLocations= cTimelapse.cTimepoint(timepoint_reg).trapLocations;
            
            xlocCELL=num2cell(xloc);
            ylocCELL = num2cell(yloc);
            
            [cTimelapse.cTimepoint(timepoint).trapLocations(:).xcenter]=deal(xlocCELL{:});
            [cTimelapse.cTimepoint(timepoint).trapLocations(:).ycenter]=deal(ylocCELL{:});
        end
        
        if ~preserve_trap_info || ~isfield(cTimelapse.cTimepoint(timepoint),'trapInfo') || isempty(cTimelapse.cTimepoint(timepoint).trapInfo)
            cTimelapse.cTimepoint(timepoint).trapInfo = trapInfo_struct;
        end
        
        % If the drift is very large then eventually the cross correlation
        % will get the wrong answer when it 'pings back' a whole trap
        % width. 
        % similarly, if the image is changing slowly over time then the
        % cross correlation will eventually become inaccurate for finding
        % drift. To prevent this, the image is replaced with the current
        % image if either image_reset_time timpoints have passed or the drift is larger
        % than half a trap width.
        if rem(timepoint-timepoints(1) +1,image_reset_time)==0 || abs(accum_row)>cTimelapse.cTrapSize.bb_height || abs(accum_col)>cTimelapse.cTrapSize.bb_width
            reg_im_fft=new_im_fft2;
            timepoint_reg=timepoint;
            accum_col = 0;
            accum_row = 0;
        end
    end
    
else
    [cTimelapse.cTimepoint(timepoints).trapLocations] = deal(struct('xcenter',0,'ycenter',0));
    if ~preserve_trap_info
        [cTimelapse.cTimepoint(timepoints).trapInfo] = deal(trapInfo_struct);
    else
        % if preserving trapInfo, only fill in those timepoints which have
        % no trapInfo.
        to_fill = cellfun('isempty',{cTimelapse.cTimepoint(timepoints).trapInfo});
        % extend if adding extra timepoints.
        to_fill(end+1:length(timepoints)) = true;
        [cTimelapse.cTimepoint(timepoints(to_fill)).trapInfo] = deal(trapInfo_struct);
    end
end

end

function [row_shift,col_shift] = dftregistration(buf1ft,buf2ft)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient image registration by crosscorrelation. Derived from code
% intended to perform sub-pixel image registration, but edited to keep only
% the whole-pixel shift part we use here (Julian Pietsch, 30 May 2022).
%
% Original code: Manuel Guizar - Dec 13, 2007
% Portions of the original code were taken from code written by Ann M. 
% Kowalczyk and James R. Fienup.
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458
% (1990).
% Citation for the original algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
% "Efficient subpixel image registration algorithms," Opt. Lett. 33,
% 156-158 (2008).
% Inputs
% buf1ft    Fourier transform of reference image,
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register,
%           DC in (1,1) [DO NOT FFTSHIFT]
% Outputs
% row_shift col_shift   Pixel shifts between images

% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
[m,n]=size(buf1ft);
CC = ifft2(buf1ft.*conj(buf2ft));
[max1,loc1] = max(CC);
[~,loc2] = max(max1);
rloc=loc1(loc2);
cloc=loc2;
md2 = fix(m/2);
nd2 = fix(n/2);
if rloc > md2
    row_shift = rloc - m - 1;
else
    row_shift = rloc - 1;
end

if cloc > nd2
    col_shift = cloc - n - 1;
else
    col_shift = cloc - 1;
end
end
