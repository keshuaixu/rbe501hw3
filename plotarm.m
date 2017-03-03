function T = plotarm( q1,q2,q3,q4,q5,q6 )

theta = deg2rad([q1;q2;q3;q4;q5;q6]);
dh = [theta(1), 290, 0, -pi/2;
    theta(2) - pi/2, 0, 270, 0;
    theta(3), 0, 70, -pi/2;
    theta(4), 302, 0, pi/2;
    theta(5), 0, 0, -pi/2;
    theta(6) - pi, 72, 0, 0];
T_intermediate = zeros(4,4,6);
for joint=1:6
    T_intermediate(:,:,joint) = dh2mat(dh(joint,1),dh(joint,2),dh(joint,3),dh(joint,4));
end

T = eye(4);
joint_cartesian = zeros([4 6]);

for joint=1:6
    T = T * T_intermediate(:,:,joint);
    joint_cartesian(:,joint) = T * [0;0;0;1];
end

joint_cartesian = horzcat([0;0;0;1], joint_cartesian)

plot3(joint_cartesian(1,:),joint_cartesian(2,:),joint_cartesian(3,:),'-o');
axis equal
view([1,-1,1])





end

