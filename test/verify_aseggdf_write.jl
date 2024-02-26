using HiQGA
vall=[[1, [2., 3.]], [2, [3., 4.]]]
sfmt = ["%3d", "%4.1f"]
outfile = "myfile"
channel_names = [["FID", "dBz/dt"], ["","pV/Am4"],["FID", "dBz/dt"]]
transD_GP.CommonToAll.writeaseggdf(vall, sfmt, outfile, channel_names)