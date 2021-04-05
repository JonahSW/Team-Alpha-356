% encke1.m        May 5, 2008

% geocentric orbit propagation using Encke's method

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global mu req j2 jdate0 tsaved rkcoef

global dtr rtd atr usaved ytbi ytrue

% astrodynamic and utility constants

om_constants;

atr = dtr / 3600;        % arc seconds to radians

% begin simulation

clc; home;

fprintf('\n                   program encke1\n');

fprintf('\n< geocentric orbit propagation using Encke''s method >\n');

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

    % request initial orbital elements

    fprintf('\ninitial orbital elements\n\n');

    oev1 = getoe([1;1;1;1;1;1]);
else
    % data file

    fprintf('\n\nplease input the simulation data file name');
    fprintf('\n(be sure to include the file name extension)\n');

    filename = input('? ', 's');

    [fid, oev1] = readoe1(filename);
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

h = 500;

% propagate orbit

fprintf('\n\n  working ...\n');

yf = rkencke ('eeqm1', neq, ti, tf, h, tetol, y);

% calculate and print final conditions

[cdstr0, utstr0] = jd2str(jdate0);

[cdstrf, utstrf] = jd2str(jdate0 + ndays);

fprintf('\n                   program encke1\n');

fprintf('\n< geocentric orbit propagation using Encke''s method >\n');

fprintf('\ninitial calendar date       ');

disp(cdstr0);

fprintf('initial universal time      ');

disp(utstr0);

fprintf('\nfinal calendar date         ');

disp(cdstrf);

fprintf('final universal time        ');

disp(utstrf);

for i = 1:1:3
    
    rf(i) = ytrue(i);
    
    vf(i) = ytrue(i + 3);
    
end

oevf = eci2orb1(mu, rf, vf);

oeprint1(mu, oevf);





