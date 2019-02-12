function gender_gap(MFNAME,WFNAME,OUTNAME)

%GENDER_GAP calculate gendergap from output files produced by MAIN.
%
%   See MAIN help for details.

% .. ingest
Tm = readtable(MFNAME);
Tf = readtable(WFNAME);

% .. get dates of each male/female file
Tm.d = datenum(Tm.x);
Tf.d = datenum(Tf.x);

% .. create standard date sequence, interp to get common observations
ds = linspace( max(Tm.d(1),Tf.d(1)), min(Tm.d(end),Tf.d(end)) );
m = interp1(Tm.d,Tm.y,ds);
f = interp1(Tf.d,Tf.y,ds);

% .. calculate gendergap, output
gender_gap = round(100*((f-m)./m)',2);
dates = datestr(ds,29);
out = table(dates,gender_gap);
writetable(out, OUTNAME)
