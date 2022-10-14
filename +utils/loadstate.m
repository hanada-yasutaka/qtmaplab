function state = loadstate(path)

   [dtype, dim, domain, hbar, basis, eval] = utils.readheader(path);
   data = utils.readdata(path, [dim, 4]);
   sysinfo = SystemInfo(dim, domain);
   vec = data(:,3) + 1i*data(:,4);
   state = FundamentalState(sysinfo, basis, vec, eval);
end

