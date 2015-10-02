
function out = sub_conv1(fin,flt,ds)

out = ifft(fin.*flt);

if ds>1
    out = out(1:ds:end)*sqrt(ds);
end

end
