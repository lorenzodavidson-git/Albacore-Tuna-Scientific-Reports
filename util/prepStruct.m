function [outputstruct] = prepStruct(sst,mld,fish,prod,type)
    
    a.lon = sst.lon;
    a.lat = sst.lat;
    a.mask = sst.mask;
    a.land = sst.land;

    mld_copy = mld;

    if type == "constant"
        a.seas = sst.constant;
        a.mld = mld_copy.constant;
        
        mld_copy.seas = mld_copy.constant; % To work in function
    elseif type == "year"
        a.seas = sst.year;
        a.mld = mld_copy.year;
        
        mld_copy.seas = mld_copy.year; % To work in function
    elseif type == "seas"
        a.seas = sst.seas;
        a.mld = mld.seas;
    end

    a = FISH_ComputeMaskGradient(a);
    a = FISH_CCS_LME_Region(a, prod);
    a = FISH_ProcessMLD(a,mld_copy);
    a = FISH_ComputeProbability_Cont(a, fish);

    outputstruct = a;