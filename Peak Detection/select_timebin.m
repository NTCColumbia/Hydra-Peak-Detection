function [time_bin] = select_timebin(hf)

time_bin = [];
count = 0;

% title('select time window, right click to stop');

while true
    if count==0
        if strcmp(get(hf,'currentkey'),'e')
           break;
        end
        current_bin = zeros(2,1);
        [coords,~,cc] = ginput(1);
        if cc==3
            break;
        end
        current_bin(1) = coords(1);
        count = count+1;
    elseif count==1
        if strcmp(get(hf,'currentkey'),'e')
           break;
        end
        [coords,~,cc] = ginput(1);
        if cc==3
            break;
        end
        current_bin(2) = coords(1);
        time_bin(end+1,:) = current_bin;
        count = 0;
    end
end

time_bin = round(time_bin);
    
end