% fill solution polygons first so they end up on the bottom
figure();
solution
[pts, shppts] = draw_users('../data/user_points');
shapedraw('../data/sites1', 'b', 'MarkerSize', 20);
shapedraw('../data/sites2', 'r', 'MarkerSize', 20);
axis([0 1920 0 1080]);
