function tofill = fillstructwithdefault(tofill, default)
% Fill a struct with default values, leaving existing fields unchanged.
%
% Inputs:
%   tofill  - Struct to be filled (may already contain some fields).
%   default - Struct containing default field-value pairs.
%
% Output:
%   tofill  - Updated struct with all missing fields from 'default' filled in.
%
    fields1 = fieldnames(tofill);
    fields2 = fieldnames(default);

    missingIdx = find(~ismember(fields2, fields1));

    for i = 1:length(missingIdx)
        tofill.(fields2{missingIdx(i)}) = default.(fields2{missingIdx(i)});
    end

end
