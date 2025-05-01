function [friction_type, ffr] = computeShellFriction(data, ratios, contact_forces, mu, dt, vel_tol)
    
    % Extract position data
    x11 = data(1:3);
    x12 = data(4:6);
    x13 = data(7:9);
    x21 = data(10:12);
    x22 = data(13:15);
    x23 = data(16:18);
    x11_0 = data(19:21);
    x12_0 = data(22:24);
    x13_0 = data(25:27);
    x21_0 = data(28:30);
    x22_0 = data(31:33);
    x23_0 = data(34:36);

    % Extract force data
    fc11 = contact_forces(1:3);
    fc12 = contact_forces(4:6);
    fc13 = contact_forces(7:9);
    fc21 = contact_forces(10:12);
    fc22 = contact_forces(13:15);
    fc23 = contact_forces(16:18);

    K2 = 15/vel_tol;
    %%
    % Compute norms
    fc11n = norm(fc11);
    fc12n = norm(fc12);
    fc13n = norm(fc13);
    fc21n = norm(fc21);
    fc22n = norm(fc22);
    fc23n = norm(fc23);

    % ratios
    r1 = ratios(1); s1 = ratios(2); t1 = ratios(3);
    r2 = ratios(4); s2 = ratios(5); t2 = ratios(6);
%% Compute relative velocities
    % tri1
    v11 = (x11 - x11_0)./dt;
    v12 = (x12 - x12_0)./dt;
    v13 = (x13 - x13_0)./dt;
    % tri2
    v21 = (x21 - x21_0)./dt;
    v22 = (x22 - x22_0)./dt;
    v23 = (x23 - x23_0)./dt;

    v1 = r1*v11 + s1*v12 + t1*v13;
    v2 = r2*v21 + s2*v22 + t2*v23;
    v_rel = v1 - v2;

    % direction of contact
    p1 = x11*r1+x12*s1+x13*t1; % point of contact on triangle 1
    p2 = x21*r2+x22*s2+x23*t2; % point of contact on triangle 2
    contact_dir = p1 - p2;
    assert(norm(contact_dir)>0, "contact distance is 0, can't compute direction of contact");
    contact_norm = contact_dir/norm(contact_dir);

    % Tangential relative velocity
    tv_rel = v_rel - (dot(v_rel, contact_norm).* contact_norm);
    tv_rel_n = norm(tv_rel);

    %%
    % Initialize output vector
    ffr = zeros(1, 18);

    % gamma: takes different value according to type of friction: 
    % zero, sliding,slipping 
    if(tv_rel_n==0)
        friction_type = "ZeroVel";
        return
    else
        if (tv_rel_n>vel_tol)
            friction_type = "Sliding";
            gamma = 1.0;
        else 
            friction_type = "Sticking";
            gamma = 2 / (1 + exp(-K2*tv_rel_n)) - 1;
        end
        
    end
    
    tv_rel_u = tv_rel / tv_rel_n;
        
    % Compute friction forces
    ffr_val = (mu*gamma).*tv_rel_u; 
    assert(~anynan(ffr_val),'IMC friction force is not real (NaN).');


    % Assign friction forces to output matrix
    ffr(1:3) = ffr_val*fc11n;
    ffr(4:6) = ffr_val*fc12n;
    ffr(7:9) = ffr_val*fc13n;
    ffr(10:12) = -ffr_val*fc21n;
    ffr(13:15) = -ffr_val*fc22n;
    ffr(16:18) = -ffr_val*fc23n;
end
