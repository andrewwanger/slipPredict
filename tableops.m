%% Load stuff
C = load('slip_data_split');
%make table
cycles = table;
names = cell(1,1426);
speed = cell(1, 1426);
nameSpeed = cell(1, 1426);
common = cell(1, 1426);
metadata = table;
for k = 1:1426
    temp = C.slip_data_split{1, k};
    spd = num2str(temp.speed) + " ";
    cyclenum = num2str(temp.cycle_id) + " ";
    newname = strcat(temp.material, spd, cyclenum);
    names{k} = temp.material;
    speed{k} = num2str(temp.speed);
    nameSpeed{k} = strcat(temp.material, num2str(temp.speed));
    common{k} = newname;
    cycles = addvars(cycles, num2cell(temp.data(1:400)),'NewVariableNames', newname);
end
varNames = cycles.Properties.VariableNames;
metadata = [metadata; names];
metadata = [metadata; speed];
metadata = [metadata; nameSpeed];
metadata.Properties.VariableNames = varNames;
combined = [cycles; metadata];
cyclesfinal = rows2vars(combined);
writetable(cyclesfinal, "cycles.csv")