function [] = rect3D(dx,dy,dz,x_0,y_0,z_0)
xarr = x_0*ones(5,6); %initialize coordinates of vertices
yarr = y_0*ones(5,6);
zarr = z_0*ones(5,6);

xarr(:,1) = xarr(:,1) + [dx; 0; 0; dx; dx]; %Parallel side to x-axis
xarr(:,2) = xarr(:,2) + dx*ones(5,1); %"Front side" in x-perspective
xarr(:,3) = xarr(:,3) + [0; 0; dx; dx; 0]; %cap
xarr(:,4) = xarr(:,1); %Opposite to side 1
xarr(:,5) = xarr(:,5); %"Back side" in x-perspective
xarr(:,6) = xarr(:,3); %base

yarr(:,1) = yarr(:,1); %"Back side" in y-perspective
yarr(:,2) =  yarr(:,2) + [0; dy; dy; 0; 0]; %Parallel side to y-axis
yarr(:,3) = yarr(:,3) + [0; dy; dy; 0; 0]; %cap
yarr(:,4) = yarr(:,1 ) + dy*ones(5,1); %"Front side" in x-perspective
yarr(:,5) = yarr(:,2); %Opposite to side 2
yarr(:,6) = yarr(:,3); %base

zarr(:,1) = zarr(:,1) + [dz; dz; 0; 0; dz]; 
zarr(:,2) = zarr(:,1);
zarr(:,3) = zarr(:,5) + dz*ones(5,1); %cap
zarr(:,4) = zarr(:,1);
zarr(:,5) = zarr(:,1); %Sides identical in z-axis

c = 2*zeros(5,1); %color set

patch(xarr, yarr, zarr, c)
end