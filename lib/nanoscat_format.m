function coeff = nanoscat_format (C, orders)
t = {};
for o = 1:length (orders);
    m = orders(o);
    
    Out = C{m};
    
    maxlen = length (Out{1});
    nlambdas = numel(Out);
    
    c = zeros(nlambdas, maxlen);
    
    for i = 1:nlambdas
        itp = interpft(Out{i}, maxlen);
        c(i, :) = itp;
    end
    t{m} = c;
end
    ct=cellfun(@transpose, t, 'uniformOutput', false);
    coeff = [ct{:}].';
end
