function [thresh] = select_thresh(hf)
% select upper and lower threshold manually on the specified figure handle


figure(hf)
[~,thresh] = ginput(2);

end