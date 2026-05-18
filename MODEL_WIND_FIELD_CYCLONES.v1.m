%__________________________________________________________________________
%
% ESTIMATION OF WIND EXPOSURE FROM HISTORICAL CYCLONE TRACKS
%
% Reconstructs past exposure to tropical cyclones in any reef location
% based on modelling the wind field around historical cyclone tracks.
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 02/2025
% Last update: 05/2026
%__________________________________________________________________________

% This code predicts the highest mean wind speed experienced at any location during a cyclone season
% based on the modelling of wind speed produced by all cyclones and the distance to cyclone centre.
% Although a location may be affected by multiple cyclones within a cyclone season (hereafter, nominal year),
% only the highest value of maximum wind speed across all cyclone events for a given season is retained.
% Therefore, the script does not report the maximum wind speed for each individual cyclone at a given location.

% Reef location data, as a table with at least the following fields
% .LAT (numeric): latitude (decimal degrees) of the location (eg, centroid of a reef polygon)
% .LON (numeric): longitude (decimal degrees) of the location

% EXAMPLE FOR THE GREAT BARRIER REEF:
MY_LOCATIONS = load('GBR_REEF_POLYGONS_2024.mat').GBR_REEFS;

% Cyclone track data need to be in the form of a table with the following fields (columns):
% (specifications for the Australian Tropical Cyclone Database)
% .NAME (categorical): name of the cyclone . Can be "Unnamed" if event was not a cyclone yet
% .TM (datetime): date/time in UTC
% .LAT, .LON (numeric): Latitude & Longitude (decimal degrees) of cyclone centre - cannot be NULL
% (sign convention: east positive, west negative & north positive, south negative)
% .CENTRAL_PRES (numeric): Central pressure of the cyclone (units: hectopascals) - can be NULL
% .ENV_PRES (numeric): Environmental pressure in which the cyclone is embedded - can be NULL
% .MN_RADIUS_MAX_WIND (numeric): mean radius (in km, from the system centre) of the maximum mean wind
% .MAX_WIND_SPD (numeric): estimated maximum mean wind in the vicinity of the cyclone centre.

% EXAMPLE FOR THE GREAT BARRIER REEF (Extracted from the Australian Tropical Cyclone Database)
TRACK_DATA = load('GBR_CYCLONE_TRACKS_2008-2026.mat').GBRMSW20082026;

% Add nominal year to the tracks
first_summer_month = 11; % November is first month of austral summer (eg, Nov 2007 -> 2008)
TRACK_DATA.NominalYear = year(TRACK_DATA.TM)+(month(TRACK_DATA.TM)>=11);

% Define maximum distance of cyclone centre to location for inclusion
max_distance = 200; % only calculate winds when distance to cyclone centre is less than 200 km

% Select wind field model:
% ModelChoice = "Boose"; % Boose et al. (2004)
ModelChoice = "Holland"; % Holland (1980) with asymmetry from McConochie et al. (1999)

OutputFileName = 'GBR_cyclones_MSW_2008-2026';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALL_YEARS = unique(TRACK_DATA.NominalYear);

Nb_loc = size(MY_LOCATIONS,1);
rho = 1.15; % air density in kg.m-3

% Preallocating the output array: we will store for every location (row) and every year (column) the highest maximum
% mean wind speed each location experienced during a cyclone season, across all cyclones that passed nearby.
% This does not record the max mean wind speed of every cyclone, only the highest maimum
PAST_MSW = nan(Nb_loc, length(ALL_YEARS));

for yr = 1:length(ALL_YEARS)

    ALL_YEARS(yr)

    I = find(TRACK_DATA.NominalYear==ALL_YEARS(yr)); % flag all cyclone positions for that nominal year
    ALL_CYCLONES_THIS_YEAR = TRACK_DATA(I,:);

    list_cyclones = unique(ALL_CYCLONES_THIS_YEAR.NAME); % list all cyclones for that year by their name.
    % Note a cyclone can be unnamed for some times - this will create duplicate tracks for the same event,
    % but not an issue here since we only keep the maximum mean wind across all events (including the "unnamed ones")

    MSW_REEF_THIS_YEAR = zeros(size(MY_LOCATIONS,1),length(list_cyclones));

    if isempty(list_cyclones)==0

        for cyc = 1:length(list_cyclones)

            list_cyclones(cyc)

            MyVar = {'NAME';'TM';'NominalYear';'LAT';'LON';'CENTRAL_PRES';'ENV_PRES';'MN_RADIUS_MAX_WIND';'MAX_WIND_SPD'};
            CURRENT_CYCLONE = ALL_CYCLONES_THIS_YEAR(ismember(ALL_CYCLONES_THIS_YEAR.NAME, list_cyclones(cyc)),MyVar);
            CURRENT_CYCLONE.VH = nan(size(CURRENT_CYCLONE,1),1); % add a column for forward motion (to be calculated)

            if size(CURRENT_CYCLONE,1)>1 % need at least 2 records to be considered

                % 1) Determine high-resolution tracking (every hour)
                CURRENT_CYCLONE_TRACK = CURRENT_CYCLONE(1,:);
                % CURRENT_CYCLONE_TRACK(1,6:end)=array2table(nan(1,5)); % initialise with NaNs
                count = 0;

                for e = 1:size(CURRENT_CYCLONE,1)-1 % for every record except the last one

                    if isnan(CURRENT_CYCLONE.MN_RADIUS_MAX_WIND(e))==1 & count==0 % ignore first records with missing radius of mximum wind
                        continue
                    else
                        count = count+1; % starts counting with the first non-NaN record
                        d = CURRENT_CYCLONE.TM(e+1)-CURRENT_CYCLONE.TM(e); % time interval with the next record

                        CURRENT_CYCLONE.VH(e) = deg2km(distance(CURRENT_CYCLONE.LAT(e), CURRENT_CYCLONE.LON(e),...
                            CURRENT_CYCLONE.LAT(e+1), CURRENT_CYCLONE.LON(e+1)))*1000/(hours(d)*3600); % hours(d) will express d in hours

                        % If radius not defined, grab the previous value
                        if isnan(CURRENT_CYCLONE.MN_RADIUS_MAX_WIND(e))==1
                            CURRENT_CYCLONE.MN_RADIUS_MAX_WIND(e) = CURRENT_CYCLONE.MN_RADIUS_MAX_WIND(e-1);
                        end

                        if hours(d) >= 2

                            lat1 = CURRENT_CYCLONE.LAT(e); lat2 = CURRENT_CYCLONE.LAT(e+1);
                            lon1 = CURRENT_CYCLONE.LON(e); lon2 = CURRENT_CYCLONE.LON(e+1);
                            [arclen,az] = distance(lat1,lon1,lat2,lon2); % arc length and azimuth of the great circle between 2 geographic positions

                            FILL = CURRENT_CYCLONE(e,:); % filling table starting with the record

                            for h = 2:hours(d)

                                FILL = [FILL ; FILL(h-1,:)]; % start by duplicating the first record
                                FILL.TM(h) = CURRENT_CYCLONE.TM(e)+hours(h-1); %increment time
                                [new_LAT, new_LON] = reckon(FILL.LAT(h-1),FILL.LON(h-1),arclen/hours(d),az); % determine virtual position every hour
                                FILL.LAT(h) = new_LAT;
                                FILL.LON(h) = new_LON;
                            end

                            CURRENT_CYCLONE_TRACK = [CURRENT_CYCLONE_TRACK ; FILL]; % add to the master table of tracking
                        end
                    end
                end

                if count~=0 % count=0 if there was no recorded value for the radius of maximum wind -> just ignore this cyclone
                    CURRENT_CYCLONE_TRACK = [CURRENT_CYCLONE_TRACK(2:end,:) ; CURRENT_CYCLONE(end,:)]; % just add the last record that was ignored while deleting the primer
                    CURRENT_CYCLONE_TRACK.VH(end) = CURRENT_CYCLONE_TRACK.VH(end-1); % extrapolate VH from previous step

                    J = find(isnan(CURRENT_CYCLONE_TRACK.MN_RADIUS_MAX_WIND)==1);
                    CURRENT_CYCLONE_TRACK.MN_RADIUS_MAX_WIND(J) = CURRENT_CYCLONE_TRACK.MN_RADIUS_MAX_WIND(J-1); % fill missing radii of maximum winds with previous ones

                    K = find(isnan(CURRENT_CYCLONE_TRACK.ENV_PRES)==1);
                    CURRENT_CYCLONE_TRACK.ENV_PRES(K)= 1010; % fill with max average pressure

                    % 2) Estimate wind for each reef along the cyclone track
                    for n = 1:length(MY_LOCATIONS.Reef_ID)

                        reef_LAT = MY_LOCATIONS.LAT(n);
                        reef_LON = MY_LOCATIONS.LON(n);

                        TRACK = CURRENT_CYCLONE_TRACK ;
                        TRACK.R = nan(size(TRACK,1),1);
                        TRACK.THETA = nan(size(TRACK,1),1);
                        TRACK.R(1) = deg2km(distance(reef_LAT, reef_LON, TRACK.LAT(1), TRACK.LON(1)));

                        % Estimate distance and angle of the reef to the cyclone's eye along track
                        for t = 2:size(CURRENT_CYCLONE_TRACK,1) % for every time position on the track
                            TRACK.R(t) = deg2km(distance(reef_LAT, reef_LON, TRACK.LAT(t), TRACK.LON(t)));
                            % Calculate azimuth (0-360 deg) of the track
                            AZ_CT = azimuth(TRACK.LAT(t-1), TRACK.LON(t-1), TRACK.LAT(t), TRACK.LON(t));
                            % Calculate azimuth (0-360 deg) of the segment between cyclone center and reef
                            AZ_CR = azimuth(TRACK.LAT(t), TRACK.LON(t), reef_LAT, reef_LON);
                            TRACK.THETA(t,1) = -(360 - AZ_CT + AZ_CR); % Counterclockwise angle between cyclone track and radial line to the reef in radians
                            % (as opposed to clockwise andgle in Northern Hemisphere - see https://github.com/hurrecon-model/HurreconR)
                        end

                        TRACK.THETA(1)=TRACK.THETA(2);  % populate THETA at t=1 (assumes same as t=2)

                        select_dist = find(TRACK.R < max_distance); % only select locations within the max distance to speed up calculations
                        SELECTION = TRACK(select_dist,4:end);

                        if isempty(select_dist)==0 % (if empty will loop to the next reef)

                            switch ModelChoice

                                case "Boose"
                                    F = 1; %scaling parameter for friction on water
                                    S = 1; % scaling parameter for asymmetry (usually set to 1)

                                    % First, calculate the Holland parameter b
                                    SELECTION.B = rho*exp(1)*(SELECTION.MAX_WIND_SPD.^2)./(100*SELECTION.ENV_PRES - 100*SELECTION.CENTRAL_PRES); % need to convert hPa in Pa
                                    A = (SELECTION.MN_RADIUS_MAX_WIND./SELECTION.R).^SELECTION.B;

                                    % Tangential wind speed in m.s-1 following Boose et al. 2004 - this is the estimated 10-min MSW on the reef
                                    SELECTION.VR = F*(SELECTION.MAX_WIND_SPD - S.*(1-sin(deg2rad(SELECTION.THETA))).*SELECTION.VH/2).*sqrt(A.*exp(1-A));

                                case "Holland"
                                    % First, calculate the Holland parameter b
                                    SELECTION.B = rho*exp(1)*(SELECTION.MAX_WIND_SPD.^2)./(100*SELECTION.ENV_PRES - 100*SELECTION.CENTRAL_PRES); % need to convert hPa in Pa
                                    A = (SELECTION.MN_RADIUS_MAX_WIND./SELECTION.R).^SELECTION.B;

                                    f = (2*7.29*1e-5)*sin(deg2rad(reef_LAT)); % Coriolis parameter at the reef

                                    % Tangential wind speed in m.s-1 following Boose et al. 2004 - this is the estimated 10-min MSW on the reef
                                    SELECTION.VR = sqrt( (SELECTION.B/rho).*A.*(100*SELECTION.ENV_PRES - 100*SELECTION.CENTRAL_PRES).*exp(-A)+(SELECTION.R*f/2).^2 ) - (SELECTION.R*f/2);

                                    % Add the effect of cyclone motion following McConochie (1999)
                                    SELECTION.alpha = 0.5*(1+cos(deg2rad(SELECTION.THETA-65)));
                                    SELECTION.VR = SELECTION.VR + SELECTION.alpha.*SELECTION.VH;
                            end

                            MSW_REEF_THIS_YEAR(n,cyc) = max(SELECTION.VR); % take the max MSW reef exposure along the track
                        end % go to next reef
                    end
                end % go to next cyclone of year yr
            end % end loop for all cyclones of year yr

            % Then only take the maximum value of MSW per reef over all the cyclones of the year
            PAST_MSW(:,yr)=max(MSW_REEF_THIS_YEAR,[],2);
        end
    end
end

%% EXPORT
save([OutputFileName '.mat'], 'PAST_MSW')
