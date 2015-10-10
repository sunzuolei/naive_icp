function err = scan_match_err(pos, scanA, scanB, thres_trans)
scanAglobal = transform_to_global(scanA, pos);
indexes=nn(scanAglobal',scanB',thres_trans);
err = sum(sum((scanAglobal(:,indexes(:,1))-scanB(:,indexes(:,2))).^2));