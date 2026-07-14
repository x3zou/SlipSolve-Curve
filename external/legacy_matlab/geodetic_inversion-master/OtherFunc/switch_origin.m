function [x_new,y_new]=switch_origin(x_old,y_old,lonc_old,latc_old,lonc_new,latc_new,ref_lon)
%% Switch the UTM origin of the data

[xo_old,yo_old]=ll2xy(lonc_old,latc_old,ref_lon);
[xo,yo]=ll2xy(lonc_new,latc_new,ref_lon);

x_new=x_old+xo_old-xo;
y_new=y_old+yo_old-yo;

end