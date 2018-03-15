# Loops through bam directory and gets read info
from __future__ import print_function
import os


dirfiles = os.listdir(".")

readsdict = {}

for file in dirfiles:
	if file.endswith(".sorted.bam.flagstat") and not file.startswith("control"):
		print(file)
		fileshortname = file.split(".")[0]
		try:
			readsdict[fileshortname]
		except KeyError:
			readsdict[fileshortname] = []
		fh = open(file, "r")
		for line in fh:
			linelist = line.strip().split(" ")
			try:
				if linelist[4] == "total":
					print(linelist)
					totalcount = linelist[0]
					if linelist[2] != "0":
						totalcount += linelist[2]
					print(totalcount)
					readsdict[fileshortname].insert(0,totalcount)  # samething a pushing to a queue
			except IndexError:
				pass # we don't care about these lines
		fh.close()

	elif file.endswith(".sorted.Q5DD.bam.flagstat") and not file.startswith("control"):
		print(file)
                fileshortname =	file.split(".")[0]
		try:
                        readsdict[fileshortname]
                except KeyError:
                        readsdict[fileshortname] = []
                
		fh = open(file, "r")
                for line in fh:
                        linelist = line.strip().split(" ")
                        try:
                            	if linelist[3] == "mapped":
                                        print(linelist)
					Q5DDreads = linelist[0]
					print(Q5DDreads)
					readsdict[fileshortname].append(Q5DDreads)
                        except IndexError:
                                pass  # we don't care all these lines
		fh.close()
	
print(readsdict)

outfh = open("readsinfo.tsv","w")

for key,value in readsdict.items():
	print("{}\t{}\t{}".format(key,value[0],value[1]))
	outfh.write("{}\t{}\t{}\n".format(key,value[0],value[1]))
outfh.close()
