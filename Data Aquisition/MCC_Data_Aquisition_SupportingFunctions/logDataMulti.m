function logDataMulti(src, event, fid)
    data = event.Data.';
    fwrite(fid, data, 'double');
end