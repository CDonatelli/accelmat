function [imu1,imu2,calib] = load_2imu(filename,calib, varargin)

opt.calib = false;
opt.clickzero = false;
opt.removeweirdblock = true;
opt.resample = true;

opt = parsevarargin(opt, varargin, 2);
%%

acc1 = h5read(filename,'/Data/Accel');
gyro1 = h5read(filename,'/Data/Gyro');
acc2 = h5read(filename,'/Data/Accel2');
gyro2 = h5read(filename,'/Data/Gyro2');
t = h5read(filename,'/Data/t');

istrigger = h5read(filename,'/Data/Trigger');
istrigger = istrigger > 0;
iszero = h5read(filename,'/Data/Zero');
iszero = iszero > 0;

acc1_units = h5readatt(filename, '/Data/Accel','Units');
acc1_units = char(acc1_units');
gyro1_units = h5readatt(filename, '/Data/Gyro','Units');
gyro1_units = char(gyro1_units');
acc1_range = h5readatt(filename, '/Data/Accel','full_range');
gyro1_range = h5readatt(filename, '/Data/Gyro','full_range');

acc2_units = h5readatt(filename, '/Data/Accel2','Units');
acc2_units = char(acc2_units');
gyro2_units = h5readatt(filename, '/Data/Gyro2','Units');
gyro2_units = char(gyro2_units');
acc2_range = h5readatt(filename, '/Data/Accel2','full_range');
gyro2_range = h5readatt(filename, '/Data/Gyro2','full_range');

t_units = h5readatt(filename, '/Data/t','Units');
t_units = char(t_units');

acc1_units = regexprep(acc1_units, '\\00', '');
gyro1_units = regexprep(gyro1_units, '\\00', '');
acc2_units = regexprep(acc2_units, '\\00', '');
gyro2_units = regexprep(gyro2_units, '\\00', '');
t_units = regexprep(t_units, '\\00', '');

%%
switch t_units
    case {'nanosec','nsec'}
        tscale = 1e9;
    case {'millisec','msec'}
        tscale = 1e3;
    otherwise
        error('Unrecognized time unit: %s\n', t_units);
end

t = t/tscale;

good = t ~= 0;
if (opt.removeweirdblock)
    t1 = t(good);
    skipind = find(diff(t1) > 1);
    if ~isempty(skipind)
        t0 = t1(skipind(end)+1);
        good = good & (t >= t0);
    end
end

t = t(good);
acc1 = acc1(good,:);
gyro1 = gyro1(good,:);
acc2 = acc2(good,:);
gyro2 = gyro2(good,:);
istrigger = istrigger(good,:);
iszero = iszero(good,:);

if any(istrigger)
    t0 = first(t,istrigger);
else
    t0 = first(t,isfinite(t));
end
t = t-t0;
%%
if opt.calib
    if opt.clickzero
        fig = figure;
        mb = msgbox('Click three zero times','Click','help','non-modal');
        plot(t,acc1);
        vertplot(t(iszero),'k--');
        
        [tzero,~] = ginputb(3);
        iszero = false(size(t));
        for i = 1:3
            [~,ind] = min(abs(t-tzero(i)));
            iszero(ind) = true;
        end
        close(mb);
        close(fig);
    end
    
    axord = inputdlg({'First zero:', 'Second zero:', 'Third zero:'}, ...
        'Order of axes',1,{'Y','-Z','X'});
    axord = regexp(axord, '([+-]?)([XYZxyz])', 'tokens', 'once');
    
    axsign = ones(3,1);
    axname = char(zeros(3,1));
    for i = 1:3
        if ~isempty(axord{i}{1}) && (axord{i}{1} == '-')
            axsign(i) = -1;
        end
        axname(i) = upper(axord{i}{2});
    end
    
    [~,ord] = sort(axname);
    axsign = axsign(ord);
    
    gax1 = acc1(iszero,:)';
    gax1 = gax1(:,ord);
    gax1 = gax1 .* repmat(axsign', [3 1]);
    
    basis1 = gramschmidt(gax1);
    %check for right handed-ness
    %the Z axis should be equal to the cross product of the X and Y axes.
    %because of small numerical issues, it's sometimes not exactly equal,
    %so we check that they're in the same direction
    assert(sum(cross(basis1(:,1),basis1(:,2)) .* basis1(:,3)) > 0.9);

    calib.chip2world1 = basis1;
    calib.world2chip1 = basis1';
end
%%
if opt.calib
    if opt.clickzero
        fig = figure;
        mb = msgbox('Click three zero times','Click','help','non-modal');
        plot(t,acc2);
        vertplot(t(iszero),'k--');
        
        [tzero,~] = ginputb(3);
        iszero = false(size(t));
        for i = 1:3
            [~,ind] = min(abs(t-tzero(i)));
            iszero(ind) = true;
        end
        close(mb);
        close(fig);
    end
    
    axord = inputdlg({'First zero:', 'Second zero:', 'Third zero:'}, ...
        'Order of axes',1,{'Y','-Z','X'});
    axord = regexp(axord, '([+-]?)([XYZxyz])', 'tokens', 'once');
    
    axsign = ones(3,1);
    axname = char(zeros(3,1));
    for i = 1:3
        if ~isempty(axord{i}{1}) && (axord{i}{1} == '-')
            axsign(i) = -1;
        end
        axname(i) = upper(axord{i}{2});
    end
    
    [~,ord] = sort(axname);
    axsign = axsign(ord);

    gax2 = acc2(iszero,:)';
    gax2 = gax2(:,ord);
    gax2 = gax2 .* repmat(axsign', [3 1]);
    
    basis2 = gramschmidt(gax2);
    %check for right handed-ness
    %the Z axis should be equal to the cross product of the X and Y axes.
    %because of small numerical issues, it's sometimes not exactly equal,
    %so we check that they're in the same direction
    assert(sum(cross(basis2(:,1),basis2(:,2)) .* basis2(:,3)) > 0.9);
    
    calib.chip2world2 = basis2;
    calib.world2chip2 = basis2';
end

%%
if opt.resample
    newsamp = diground(1./nanmean(diff(t)), 10);
    fprintf('Resampling at %d Hz\n', newsamp);
    t0 = t;
    t = (t0(1):1/newsamp:t0(end))';
    acc1 = interp1(t0,acc1, t);
    gyro1 = interp1(t0,gyro1, t);
    acc2 = interp1(t0,acc2, t);
    gyro2 = interp1(t0,gyro2, t);
end
%% 
imu1.t = t;
imu1.acc = acc1 * calib.chip2world2;
imu1.acc_units = acc1_units;
imu1.acc_range = acc1_range;
imu1.gyro = gyro1 * calib.chip2world1;
imu1.gyro_units = gyro1_units;
imu1.gyro_range = gyro1_range;

imu2.t = t;
imu2.acc = acc2 * calib.chip2world2;
imu2.acc_units = acc2_units;
imu2.acc_range = acc2_range;
imu2.gyro = gyro2 * calib.chip2world2;
imu2.gyro_units = gyro2_units;
imu2.gyro_range = gyro2_range;