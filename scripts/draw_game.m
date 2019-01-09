% fill solution polygons first so they end up on the bottom
figure();
solution
[pts, shppts] = draw_users('../data/user_points');
shapedraw('../data/sites1', 'b', 'MarkerSize', 20);
shapedraw('../data/sites2', 'r', 'MarkerSize', 20);
% draw the solution as a double-circle for player 2
scatter(spoint(:,1), spoint(:,2), [600], 'w', 'filled');
scatter(spoint(:,1), spoint(:,2), [400], 'r', 'filled');
axis([0 1920 0 1080]);
