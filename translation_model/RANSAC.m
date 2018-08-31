function T = RANSAC(f1,f2,matches)
    x1 = f1(1,matches(1,:)); x2 = f2(1,matches(2,:));
    y1 = f1(2,matches(1,:)); y2 = f2(2,matches(2,:));
    dx = x2-x1; dy = y2-y1;
    [density_x,xi] = ksdensity(dx); [density_y,yi] = ksdensity(dy);
    [~,dx_opt_ind]= max(density_x); [~,dy_opt_ind]= max(density_y);
    dx_opt = xi(dx_opt_ind); dy_opt = yi(dy_opt_ind);
    dx_opt = round(dx_opt); dy_opt = round(dy_opt);
    T = [abs(dx_opt),dy_opt];
end