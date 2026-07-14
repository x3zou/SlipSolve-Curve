function area = triangle_area(v1, v2, v3)
    % Xiaoyu Zou, 2/28/2025
    % TRIANGLE_AREA Computes the area of a triangle in 3D space
    % given its three vertices v1, v2, and v3.
    %
    % Input:
    %   v1, v2, v3 - 1x3 vectors representing the coordinates of the vertices
    %
    % Output:
    %   area - scalar value of the triangle's area

    % Compute two edge vectors
    edge1 = v2 - v1;
    edge2 = v3 - v1;

    % Compute the cross product of the two edge vectors
    cross_product = cross(edge1, edge2);

    % Compute the area using the norm of the cross product
    area = 0.5 * norm(cross_product);
end