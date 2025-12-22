function X = getSeq(Xin, layout)
if strcmpi(layout,'KxT')
    X = Xin;
elseif strcmpi(layout,'TxK')
    X = Xin';
else
    error('opt.dataLayout must be ''KxT'' or ''TxK''.');
end
end