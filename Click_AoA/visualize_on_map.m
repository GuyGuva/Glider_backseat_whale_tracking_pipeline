function [lat_dec,lon_dec]=visualize_on_map(d,clickBearing,co,lat,lon,head)

% Convert DMS to decimal degrees
lat_str=num2str(lat);  lon_str=num2str(lon);
Deg_lat=str2num(lat_str(1:2));  Deg_lon=str2num(lon_str(2:3));
Min_lat=str2num(lat_str(3:4));  Min_lon=str2num(lon_str(4:5));
Sec_lat=str2num(lat_str(6:end));  Sec_lon=str2num(lon_str(7:end));
lat_dec=Deg_lat+Min_lat/60+Sec_lat/3600;
lon_dec=-(Deg_lon+Min_lon/60+Sec_lon/3600);


    lat=lat_dec*ones(size(clickBearing));
    lon=lon_dec*ones(size(clickBearing));

    % Earth radius and line length
    R = 6371000;       % [m]
    % d = 1000;          % [m]  % 1 km
    d_arrow = 500;      % orientation arrow length [m]
    d_head   = 200;        % arrowhead side length [m]
    theta    = deg2rad(30); % arrowhead opening angle (30°)
    
    % Convert to radians
    lat1 = deg2rad(lat);
    lon1 = deg2rad(lon);
    brg  = deg2rad(clickBearing(:));  % ensure column

    brg_az  = deg2rad(head(:));          % sensor azimuth

    
    % Destination point after distance d along bearing brg
    lat2 = asin( sin(lat1).*cos(d/R) + cos(lat1).*sin(d/R).*cos(brg) );
    lon2 = lon1 + atan2( sin(brg).*sin(d/R).*cos(lat1), ...
                         cos(d/R) - sin(lat1).*sin(lat2) );
    
    % Back to degrees
    lat2 = rad2deg(lat2);
    lon2 = rad2deg(lon2);
%%

lat_tip = asin( sin(lat1).*cos(d_arrow/R) + ...
                cos(lat1).*sin(d_arrow/R).*cos(brg_az) );

lon_tip = lon1 + atan2( sin(brg_az).*sin(d_arrow/R).*cos(lat1), ...
                        cos(d_arrow/R) - sin(lat1).*sin(lat_tip) );

lat_tip_deg = rad2deg(lat_tip);
lon_tip_deg = rad2deg(lon_tip);

% Bearings (from tip) for left/right arrowhead segments
brg_left  = brg_az + pi - theta;
brg_right = brg_az + pi + theta;

% Left head point
lat_headL = asin( sin(lat_tip).*cos(d_head/R) + ...
                  cos(lat_tip).*sin(d_head/R).*cos(brg_left) );
lon_headL = lon_tip + atan2( sin(brg_left).*sin(d_head/R).*cos(lat_tip), ...
                             cos(d_head/R) - sin(lat_tip).*sin(lat_headL) );

% Right head point
lat_headR = asin( sin(lat_tip).*cos(d_head/R) + ...
                  cos(lat_tip).*sin(d_head/R).*cos(brg_right) );
lon_headR = lon_tip + atan2( sin(brg_right).*sin(d_head/R).*cos(lat_tip), ...
                             cos(d_head/R) - sin(lat_tip).*sin(lat_headR) );

lat_headL_deg = rad2deg(lat_headL);
lon_headL_deg = rad2deg(lon_headL);
lat_headR_deg = rad2deg(lat_headR);
lon_headR_deg = rad2deg(lon_headR);


%%
% figure;
subplot(1,3,3);
ax  = gca;               % Get the axes created by subplot
pos = ax.Position;       %  Save the position BEFORE deleting
delete(ax);              % Remove the normal axes

gx = geoaxes('Position', pos);   %  Use 'pos', not ax.Position
hold(gx, 'on');

% Sensor positions (points)
geoscatter(gx, lat, lon, 20, 'filled', 'MarkerFaceColor', 'k'); hold on;

% For each click: draw arrow shaft + arrowhead
for k = 1:numel(lat)
    % Shaft (sensor position → tip)
    geoplot(gx, [lat(k) lat_tip_deg(k)], ...
                [lon(k) lon_tip_deg(k)], ...
                'r-', 'LineWidth', 1.5);

    % Arrowhead: two small lines from tip to left/right head points
    geoplot(gx, [lat_tip_deg(k) lat_headL_deg(k)], ...
                [lon_tip_deg(k) lon_headL_deg(k)], ...
                'r-', 'LineWidth', 1.2);
    geoplot(gx, [lat_tip_deg(k) lat_headR_deg(k)], ...
                [lon_tip_deg(k) lon_headR_deg(k)], ...
                'r-', 'LineWidth', 1.2);
end


    for k = 1:numel(lat)
        Color=co{k};
        geoplot(gx, [lat(k) lat2(k)], [lon(k) lon2(k)],Color(1)); hold on;
    end

    title(gx, 'Sensor position and estimated source bearings');

% Sensor positions (points) – static

% geoscatter(gx, lat, lon, 20, 'filled', 'MarkerFaceColor', 'k');
% 
% % ---- Create empty line objects ONCE ----
% hShaft  = geoplot(gx, NaN, NaN, 'r-', 'LineWidth', 1.5);  % azimuth shaft
% hHeadL  = geoplot(gx, NaN, NaN, 'r-', 'LineWidth', 1.2);  % arrowhead left
% hHeadR  = geoplot(gx, NaN, NaN, 'r-', 'LineWidth', 1.2);  % arrowhead right
% hBearing = geoplot(gx, NaN, NaN, 'b-', 'LineWidth', 1.0); % 1 km source-bearing
% 
% title(gx, 'Sensor position and estimated source bearings');
% 
% % ---- Loop over clicks and UPDATE the same lines ----
% 
% for k = 1:numel(lat)
% 
%     % 1) Update arrow shaft (sensor → azimuth tip)
%     set(hShaft, 'LatitudeData',  [lat(k)      lat_tip_deg(k)], ...
%                 'LongitudeData', [lon(k)      lon_tip_deg(k)]);
% 
%     % 2) Update arrowhead (tip → left/right head points)
%     set(hHeadL, 'LatitudeData',  [lat_tip_deg(k) lat_headL_deg(k)], ...
%                 'LongitudeData', [lon_tip_deg(k) lon_headL_deg(k)]);
% 
%     set(hHeadR, 'LatitudeData',  [lat_tip_deg(k) lat_headR_deg(k)], ...
%                 'LongitudeData', [lon_tip_deg(k) lon_headR_deg(k)]);
% 
% 
%     % 3) Update 1 km source-bearing line (sensor → source direction)
%     set(hBearing, 'LatitudeData',  [lat(k) lat2(k)], ...
%                   'LongitudeData', [lon(k) lon2(k)]);
% 
% 
% end
%     drawnow;           % update the geoaxes
%     pause(0.05);       % slow down the animation


end