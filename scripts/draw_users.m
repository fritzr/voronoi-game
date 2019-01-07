function [points, shppoints, rings, shprings] ...
    = draw_users (filename, draw_index, axes)
% -- draw_users(FILE[, DRAW_INDEX[, AXES]])
% -- POINTS, RINGS, SHP_POINTS, SHP_RINGS = draw_users(...)
%     Read, draw, and return the points and rings read from the ESRI
%     shapefiles named <FILE>_points.{shp,dbf} and <FILE>_rings.{shp,dbf}.
%
%     DRAW_INDEX controls which user (and associated rings) to draw, if any:
%     'none'               -- do not plot anything
%     scalar or 1xN vector -- draw each user with matching index (1-based)
%     'all'                -- draw all users (default if omitted)
%
%     If DRAW_INDEX is not 'none' and AXES is given, use axis(AXES) for plots.
%
%     The return vectors are as follows:
%     POINTS     -- Nx2 vector of [ X Y ] center user points.
%     SHP_POINTS -- struct array of all points
%     RINGS      -- 1xM cell array: each element is an Nx2 polygon vector
%     SHP_RINGS  -- struct array of all rings
%
%     If RINGS is omitted, only POINTS will be read/drawn and no modifications
%     to the filename will be used (we will read <FILE>.{shp,dbf}).
%
%     For help on the format of SHP_POINTS/SHP_RINGS, see the shaperead
%     function. This function requires that the geometry, and shapefile
%     packages be installed.
%
%     See also: shaperead, shapedraw.

  if nargin < 2,
    draw_index = 'all';
  end
  if nargin < 3,
    axes = [0 1920 0 1080];
  end

  read_rings = nargout >= 3;

  if read_rings,
    point_file = [filename '_points'];
    rings_file = [filename '_rings'];
  else
    point_file = filename ;
    rings_file = '';
  end

  shppoints = shaperead(point_file);
  points = [ [shppoints.X]' [shppoints.Y]' ];

  switch draw_index
    case 'none'
      draw_index = [];
    case 'all'
      draw_index = colon(1, length(points));
  end

  if read_rings,
    shprings = shaperead(rings_file);
    nrings = length(shprings);
    rings = cell(1, nrings);
    for i = 1:nrings,
      rings{i} = [ [shprings(i).X]' [shprings(i).Y]' ];
    end
  end

  if length(draw_index) > 0,
    cmap = colormap();
    ncolors = length(cmap);
    cm_idxmap = linspace(1, ncolors, length(draw_index));

    widths = 200 .* ones(1,length(draw_index));
    scatter(points(draw_index,1), points(draw_index,2), widths, draw_index, ...
        'filled');

    if ! read_rings,
      return
    end

    % draw rings with colors matching the center point,
    for i = 1:length(rings),
      % but only if the ring's point is in draw_index
      loc = find(draw_index == (shprings(i).pointIndex + 1));
      if length(loc) > 0,
        color = cmap(cm_idxmap(loc), :);
        phnd = patch(rings{i}(:,1), rings{i}(:,2), color);
        set(phnd, 'FaceAlpha', 0, 'FaceColor', 'none', 'EdgeColor', color, ...
            'LineWidth', 2);
      end
    end

    axis(axes);
  end

end
