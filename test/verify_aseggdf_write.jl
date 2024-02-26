using HiQGA
vall=[[100, 1, [2., 3.]], [100, 2, [3., 4.]]]
sfmt = ["%4i", "%3i", "%4.1f"]
outfile = "myfile"
channel_names = [["Line", "FID", "dBzdt"], ["", "","pV/Am4"],["Line", "FID", "dBz/dt"]]
transD_GP.CommonToAll.writeasegdat(vall, sfmt, outfile, channel_names)
transD_GP.CommonToAll.writeasegdfn(vall, channel_names, sfmt, outfile)
