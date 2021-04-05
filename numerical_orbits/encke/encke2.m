% encke2.m        May 5, 2008

% heliocentric orbit propagation using Encke's method

% perturbations include Venus, Earth, Mars, Jupiter and Saturn

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global rtd dtr

global mu pmu jdate0 tsaved rkcoef

global aunit atr usaved ytbi ytrue

% astrodynamic and utility constants

om_constants;

% arc seconds to radians

atr = dtr / 3600;

% gravitational constants (km^3/sec^2)

pmu(1) = 132712438000;        % sun
pmu(2) = 324858.77;           % venus
pmu(3) = 398600.485;          % earth
pmu(4) = 42828.287;           % mars
pmu(5) = 126685807.786;       % jupiter
pmu(6) = 37924490.093;        % saturn

mu = pmu(1);

% begin simulation

clc; home;

fprintf('\n                    program encke2\n');

fprintf('\n< heliocentric orbit propagation using Encke''s method >\n');

% request simulation data

while (1)
    
    fprintf('\n\n simulation data menu\n');

    fprintf('\n   <1> user input\n');

    fprintf('\n   <2> data file\n\n');

    slct = input('? ');

    if (slct == 1 || slct == 2)
        break;
    end
    
end

if (slct == 1)
    
    % request initial calendar date and UT

    fprintf('\n\ninitial calendar date and universal time\n');

    [month, day, year] = getdate;

    [uthr, utmin, utsec] = gettime;

    while (1)
        fprintf('\nplease input the simulation period (days)\n');

        ndays = input('? ');

        if (ndays > 0)
            break;
        end
        
    end

    while (1)
        
        fprintf('\nplease input the algorithm error tolerance');
        fprintf('\n(a value between 1.0e-8 and 1.0e-10 is recommended)\n');

        tetol = input('? ');

        if (tetol > 0)
            break;
        end
        
    end

    % total simulation time (seconds)

    tf = 86400 * ndays;

    % initial julian date

    jdate0 = julian(month, day, year) ...
        + uthr / 24 + utmin / 1440 + utsec / 86400;

    % initial orbit

    fprintf('\ncalendar date and universal time of perihelion passage\n');

    [month, day, year] = getdate;

    [uthr, utmin, utsec] = gettime;

    jdpp = julian(month, day, year) ...
        + uthr / 24 + utmin / 1440 + utsec / 86400;

    while (1)
        
        fprintf('\n\nplease input the perihelion radius (AU''s)');
        fprintf('\n(perihelion radius > 0)\n');

        rper = input('? ');

        if (rper > 0)
            break;
        end
        
    end

    rper = aunit * rper;

    % request additional orbital elements

    fprintf('\norbital elements - true ecliptic and equinox of date\n\n');

    oev1 = getoe([0;1;1;1;1;0]);

    ecc = oev1(2);

    sma = rper / (1 - ecc);

    % time since perhelion passage (seconds)

    tspp = 86400 * (jdate0 - jdpp);

    % compute initial true anomaly

    manom = sqrt(mu(1) / abs(sma) ^ 3) * tspp;

    [eanom, tanom] = kepler1 (manom, ecc);

    oev1(1) = sma;
    
    oev1(6) = tanom;

else
    
    % data file

    fprintf('\n\nplease input the simulation data file name');
    fprintf('\n(be sure to include the file name extension)\n');

    filename = input('? ', 's');

    [fid, oev1, jdpp] = readoe2(filename);

    fprintf('\n\ninitial calendar date and universal time\n');

    [month, day, year] = getdate;

    [uthr, utmin, utsec] = gettime;

    while (1)
        
        fprintf('\nplease input the simulation period (days)\n');

        ndays = input('? ');

        if (ndays > 0)
            break;
        end
        
    end

    while (1)
        
        fprintf('\nplease input the algorithm error tolerance');
        fprintf('\n(a value between 1.0e-8 and 1.0e-10 is recommended)\n');

        tetol = input('? ');

        if (tetol > 0)
            break;
        end
        
    end

    % initial julian date

    jdate0 = julian(month, day, year) ...
        + uthr / 24 + utmin / 1440 + utsec / 86400;

    % total simulation time (seconds)

    tf = 86400 * ndays;

    rper = aunit * oev1(1);

    ecc = oev1(2);

    sma = rper / (1 - ecc);

    % time since perhelion passage (seconds)

    tspp = 86400 * (jdate0 - jdpp);

    % compute initial true anomaly

    manom = sqrt(mu(1) / abs(sma) ^ 3) * tspp;

    [eanom, tanom] = kepler1 (manom, ecc);

    oev1(1) = sma;
    
    oev1(6) = tanom;
    
end

% calculate initial position and velocity vectors

[rtbi, vtbi] = orb2eci(mu, oev1);

% initialize deviation (y) and reference (ytbi) vectors

for i = 1:1:6
    
    y(i) = 0;
    
end

for i = 1:1:3
    
    ytbi(i) = rtbi(i);
    
    ytbi(i + 3) = vtbi(i);
    
end

% initialize

ti = 0;

tsaved = 0;

% initialize rkf integrator

rkcoef = 1;

usaved = 0;

neq = 6;

h = 50 * 86400;

% propagate orbit

yf = rkencke ('eeqm2', neq, ti, tf, h, tetol, y);

% calculate and print final conditions

[cdstr1, utstr1] = jd2str(jdate0 + ndays);

fprintf('\n                    program encke2\n');

fprintf('\n< heliocentric orbit propagation using Encke''s method >\n');

fprintf('\n\nfinal conditions\n');

fprintf('\ncalendar date      ');

disp(cdstr1);

fprintf('\nuniversal time     ');

disp(utstr1);

for i = 1:1:3
    
    rf(i) = ytrue(i);
    
    vf(i) = ytrue(i + 3);
    
end

oevf = eci2orb1(mu, rf, vf);

oeprint2(mu, oevf);

rper = oevf(1) * (1 - oevf(2));

fprintf('\n\nperihelion radius   %12.10e AU \n\n', rper/aunit);




