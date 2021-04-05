function ydot = eeqm2 (time, y)

% first order equations of heliocentric orbital motion

% required by encke2.m

% input

%  time = simulation time (seconds)
%  y    = state vector

% output

%  ydot = integration vector

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global pmu jdate0 tsaved ytbi

jdate = jdate0 + time /86400;

tau = time - tsaved;

if (tau ~= 0)
    
    % compute current two-body orbit

    ytbf = twobody4 (pmu(1), tau, ytbi);

    % compute current true orbit

    for i = 1:1:6
        
        ytrue(i) = ytbf(i) + y(i);
        
    end

    % calculate "true" and "osculating" heliocentric position

    rs2sc = norm(ytrue(1:3));

    rrs2sc = pmu(1) / rs2sc ^ 3;

    rtb = norm(ytbf(1:3));

    rrtb = pmu(1) / rtb ^ 3;

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

% calculate planetary state vectors

% Venus

rvenus = venus(jdate);

rpm(1) = norm(rvenus);

for i = 1:1:3
    
    rp(1, i) = rvenus(i);
    
end

% Earth

rearth = earth(jdate);

rpm(2) = norm(rearth);

for i = 1:1:3
    
    rp(2, i) = rearth(i);
    
end

% Mars

rmars = mars(jdate);

rpm(3) = norm(rmars);

for i = 1:1:3
    
    rp(3, i) = rmars(i);
    
end

% Jupiter

rjupiter = jupiter(jdate);

rpm(4) = norm(rjupiter);

for i = 1:1:3
    
    rp(4, i) = rjupiter(i);
    
end

% Saturn

rsaturn = saturn(jdate);

rpm(5) = norm(rsaturn);

for i = 1:1:3
    
    rp(5, i) = rsaturn(i);
    
end

for i = 1:1:5
    
    rrpm(i) = pmu(i + 1) / (rpm(i) ^ 3);
    
end

% compute planetocentric position vectors of spacecraft

for i = 1:1:5
    
    rp2sc(i, 1) = ytrue(1) - rp(i, 1);
    
    rp2sc(i, 2) = ytrue(2) - rp(i, 2);
    
    rp2sc(i, 3) = ytrue(3) - rp(i, 3);
    
end

% compute planetocentric distances of spacecraft

for i = 1:1:5
    
    rp2scm(i) = sqrt(rp2sc(i, 1) ^ 2 + rp2sc(i, 2) ^ 2 + rp2sc(i, 3) ^ 2);

    rrp2scm(i) = pmu(i + 1) / (rp2scm(i) ^ 3);
    
end

% compute planetary perturbations

for j = 1:1:3

    accp(j) = 0;

    for i = 1:1:5
        
        accp(j) = accp(j) - rp2sc(i, j) * rrp2scm(i) - rp(i, j) * rrpm(i);
        
    end
end

% compute integration vector

ydot(1) = y(4);

ydot(2) = y(5);

ydot(3) = y(6);

if (tau ~= 0)
    
    ydot(4) = -rrtb * (f * ytrue(1) + y(1)) + accp(1);

    ydot(5) = -rrtb * (f * ytrue(2) + y(2)) + accp(2);

    ydot(6) = -rrtb * (f * ytrue(3) + y(3)) + accp(3);
    
else
    
    % time = 0 - only planetary perturbations

    ydot(4) = accp(1);

    ydot(5) = accp(2);

    ydot(6) = accp(3);
    
end
