function [sunrise,sunset] = sunTimes(lat,lon,dates,tz,useStandardTime)

% Code from: https://www.mathworks.com/matlabcentral/fileexchange/55509-sunrise-sunset

% useStandardTime = false;

yearStart = dateshift(dates(1,1),"start","year");               % midnight beginning of year
yearEnd = dateshift(dates(1,1),"start","year","next");          % midnight end of year
noon = dateshift(dates,"start","day") + hours(12);           % noon on each date of interest

% calculate Equation of Time
B = 360*(day(noon,'dayofyear')-81)/365;
eot = minutes(9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B));

% calculate Hour Angle (omega)
yearFrac = (noon - yearStart) ./ (yearEnd - yearStart);    % year fraction at noon each day
gamma = 2*pi * yearFrac;                                   % year fraction in radians
delta = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 0.006758*cos(2*gamma) ...
    + 0.000907*sin(2*gamma) - 0.002697*cos(3*gamma) + 0.00148*sin(3*gamma);
omega = acosd((cosd(90.833)./(cosd(lat).*cos(delta))) - tand(lat).*tan(delta));

% calculate sunrise & sunset
sunrise = noon - minutes(4*(lon + omega)) - eot;
sunset  = noon - minutes(4*(lon - omega)) - eot;

% tz = locationToTimeZone(lat,lon,useStandardTime)
% sunrise.TimeZone = tz;
% sunset.TimeZone = tz;

% adjust for time zone
if ~isempty(tz) && tz~=0
    sunrise = sunrise + hours(tz);
    sunset = sunset + hours(tz);
end

end


function tz = locationToTimeZone(lat,lon,useStandardTime)
load timeZones.mat timeZones

lat = single(lat);
lon = single(lon);
for i = 1:size(timeZones,1)
    if inpolygon(lon, lat, timeZones(i).Lon,timeZones(i).Lat)
        if useStandardTime
            tz = timeZones(i).FixedID;
        else
            tz = timeZones(i).ID;
        end
        return
    end
end

tz = compose("Etc/GMT%+d",-floor((lon+7.5)/15)); % opposite sign conventions

end
