using HiQGA
vall=[[1, [2., 3.]], [2, [3., 4.]]]
sfmt = ["%3i", "%4.1f"]
outfile = "myfile"
channel_names = [["FID", "dBzdt"], ["","pV/Am4"],["FID", "dBz/dt"]]
transD_GP.CommonToAll.writeaseggdf(vall, sfmt, outfile, channel_names)
transD_GP.CommonToAll.writedfn(vall, channel_names, sfmt, outfile)
