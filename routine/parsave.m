function [ ] = parsave( name, result )
% save the results in parfor
    save(name, 'result');
end

