function normTime = normTimeofD(times,sunrise,sunset)
% From A. Solsona-Berga, ~2019
% Given a set of detections and diel informaiton specifying night time,
% renormalize detections to represent a 12 hour day/night period by
% linear interpolation.
%
% Assumptions:
% Both detections and night are sorted by timestamp and
%   converted to local time (or in UTC with a provided UTCoffset)
%   so that night fall is after sunrise each day.
% There are no detections outside of the night intervals except for
%   the day before and after the first and last night respectively.

N = length(times);
normTime = nan(N,1); % initialize

% 1 where there is sunrise and 0 where there is sunset
sunsetFirst = sunrise(1)>sunset(1);
if sunsetFirst
    sunset = [sunset(2:end);sunset(end)];
end

sunrise_next = [sunrise(2:end);sunrise(end)];

totalDay = minutes(sunset-sunrise);
totalNight = minutes(sunrise_next-sunset);
%process each time: and normalize by the maximum minutes between rise/set

if times(1) < sunrise
    % No diel information for this date
    error('%s earlier than %s. Add previous diel information', ...
        datestr(times(1)), datestr(sunrise));
end
if times(end) > sunset
    % No diel information for this date
    error('%s earlier than %s. Add posterior diel information', ...
        datestr(times(end)), datestr(sunset));
end

for i = 1:N
    event = times(i);
    if sum(event >= sunrise & event < sunset)
        % day time, convert to normalized day
        idxDay = find(event >= sunrise & event < sunset);
        normTime(i) = - minutes(sunset(idxDay)-event)/totalDay(idxDay);
    elseif sum(event >= sunset & event < sunrise_next)
        % night time, convert to normalized night
        idxNight = find(event >= sunset & event < sunrise_next);
        normTime(i) = minutes(event-sunset(idxNight))/totalNight(idxNight);
    else
        disp('Warning: Detection time outside range of sunrise/sunset intervals provided, moving on to next detection\n')
    end
end
end