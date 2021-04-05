function ydot = eeqm1 (time, y)

% first order equations of orbital motion

% required by encke1.m

% input

%  time = simulation time (seconds)
%  y    = state vector

% output

%  ydot = integration vector

% global

%  ytbi = "two body" state vector

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu req j2 tsaved ytbi

tau = time - tsaved;

if (tau ~= 0)
    
    % compute current two-body orbit

    ytbf = twobody4 (mu, tau, ytbi);

    % compute current true orbit

    for i = 1:1:6
        
        ytrue(i) = ytbf(i) + y(i);
        
    end

    % calculate "true" and "osculating" position

    rs2sc = sqrt(ytrue(1) * ytrue(1) + ytrue(2) * ytrue(2) + ytrue(3) * ytrue(3));

    rrs2sc = mu / rs2sc ^ 3;

    rtb = sqrt(ytbf(1) * ytbf(1) + ytbf(2) * ytbf(2) + ytbf(3) * ytbf(3));

    rrtb = mu / rtb ^ 3;

    % compute f(q)

    q = 0;

    for i = 1:1:3
        
        q = q + y(i) * (y(i) - 2 * ytrue(i));
        
    end

    q = q / (rs2sc * rs2sc);

    q1 = (q + 1) ^ 1.5;

    f = q * (3 + 3 * q + q * q) / (1 + q1);
    
else
    
    ytrue = ytbi;
    
end

% compute j2 perturbations

r2 = ytrue(1) * ytrue(1) + ytrue(2) * ytrue(2) + ytrue(3) * ytrue(3);

r1 = sqrt(r2);

r3 = r2 * r1;

r5 = r2 * r3;

d1 = -1.5 * j2 * req * req * mu / r5;

d2 = 1 - 5 * ytrue(3) * ytrue(3) / r2;

% j2-only gravity acceleration vector

agrav(1) = ytrue(1) * d1 * d2;

agrav(2) = ytrue(2) * d1 * d2;

agrav(3) = ytrue(3) * d1 * (d2 + 2);

% compute integration vector

ydot(1) = y(4);

ydot(2) = y(5);

ydot(3) = y(6);

if (tau ~= 0)
    
    ydot(4) = -rrtb * (f * ytrue(1) + y(1)) + agrav(1);

    ydot(5) = -rrtb * (f * ytrue(2) + y(2)) + agrav(2);

    ydot(6) = -rrtb * (f * ytrue(3) + y(3)) + agrav(3);
    
else
    
    % time = 0 ==> perturbations only

    ydot(4) = agrav(1);

    ydot(5) = agrav(2);

    ydot(6) = agrav(3);
    
end
