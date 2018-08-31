function [M,T] = RANSAC(matches,f_obj,f_sce,delta_th)
    M = zeros(2,2);
    T = zeros(2,1);
    n_iter = 1000; %number of iterations
    n_match = size(matches,2);
    n_points = 3;
    n_corrisp_max = 0; %number of corrispondence maximum
    
    for n=1:n_iter
        n_corrisp = 0; %number of corrispondence
        obj_ind = round(rand(1,n_points)*(n_match-1)+1);
        sce_ind = matches(2,obj_ind);
        obj_points = f_obj(1:2,obj_ind);
        sce_points = f_sce(1:2,sce_ind);
        b = zeros(2*n_points,1);
        for i=1:n_points
            b(2*(i-1)+1) = sce_points(1,i);
            b(2*(i-1)+2) = sce_points(2,i);
        end
        
        A = zeros(2*n_points,6);
        for i=1:n_points
            A(2*(i-1)+1,:) = [obj_points(1,i),obj_points(2,i),0,0,1,0];
            A(2*(i-1)+2,:) = [0,0,obj_points(1,i),obj_points(2,i),0,1];
        end
        x_temp = inv(A'*A)*A'*b;
        M_temp = [x_temp(1),x_temp(2);x_temp(3),x_temp(4)];
        T_temp = [x_temp(5);x_temp(6)];
        
        for i=1:n_match
            for j=1:n_points
                if i~=sce_ind(j);
                    obj_point = f_obj(1:2,i);
                    sce_point = M_temp*obj_point + T_temp;
                    du = sce_point(1)-f_sce(1,matches(2,i));
                    dv = sce_point(2)-f_sce(2,matches(2,i));
                    delta = abs(du)+abs(dv);
                    %disp(delta)

                    if delta<delta_th
                        n_corrisp = n_corrisp+1;
                        if n_corrisp > n_corrisp_max
                            n_corrisp_max = n_corrisp;
                            M = M_temp;
                            T = T_temp;
                        end
                    end 
                end
            end
        end
        
    end
end