#!/usr/bin/env python
# MinPath (current version: 1.21)
# developed by Yuzhen Ye (yye@indiana.edu)
# Indiana University, Bloomington

# version 1.21 fixed a bug that messed up with fig input (released on Jan 11, 2018)
# what's new in version 1.2 (released on Oct 19, 2010)
#  MinPath1.2 works on any pathway system
#  what you need: a pathway-function mapping file (e.g, data/ec2path), and your input file

# what's new in version 1.1 (released on July 31, 2010) 
#  * numpy independent
#  * add detailed output of ko assignments to each pathway

import sys
import os
import re
import operator
import math

###path0 patch added by Fernando Puente-Sanchez. 08-VII-2018
path0 = os.path.dirname(os.path.realpath(sys.argv[0]))
minpath = path0

#minpath = os.environ.get('MinPath')
#path0 = "/home/fer/opt/MinPath"
#if minpath or os.path.exists(path0):
#        if os.path.exists(path0):
#                minpath = path0
#else:
#        sys.exit("Environment variable MinPath not set")

keggPath0, seedPath0, mapPath0, glpsol0 = minpath + "/../data", minpath + "/../data", minpath + "/../data", minpath + "/glpsol"

#################



def intmatrix(dim1, dim2):
	mat = []
	for i in range(dim1):
		tmp = [0] * dim2		
		mat.append(tmp)

	return mat
	
class MinPath:

	def __init__(self, whichdb = "KEGG", pathdir="", getgene = False, givenspe = "", mapfile = ""):
		self.whichDB = whichdb
		self.speID = givenspe
		global keggPath0
		global seedPath0
		global mapPath0
		if pathdir == "":
			if whichdb == "KEGG": 
				self.dataDir = keggPath0
			elif whichdb == "SEED":
				self.dataDir = seedPath0
			else:
				self.dataDir = mapPath0
		else:
			self.dataDir = pathdir 
		self.famTot = 0
		self.pathwayTot = 0
		self.famList = []
		self.famID = []
		self.famName = []
		self.famCount = []
		self.fam2Path = []
		self.pathList = []
		self.pathName = []
		self.path2Fam = []
		self.path2FamUni = []
		self.orgList = []
		self.orgGene2Fam = []
		self.orgGeneList = []
		self.famMapped = []
		self.pathMapped = []
	
		if whichdb == "SEED":
			print "now get SEED"
			fig2ssfile = self.dataDir + "/figfam_subsystem.dat"
			self.ReadFigSubsytem(fig2ssfile)
		elif whichdb == "KEGG":
			pathfile = self.dataDir + "/map_title.tab"	
			self.ReadKEGGPath(pathfile)

			kofile = self.dataDir + "/ko" 
			self.ReadKO(kofile, getgene, givenspe) #KO: KEGG family
		else:
			if os.path.exists(mapfile):
				print "mapfile", mapfile
				self.ReadAnyMap(mapfile)
			elif os.path.exists(self.dataDir + "/" + os.path.basename(mapfile)):
				print "file: ", self.dataDir + "/" + os.path.basename(mapfile)
				self.ReadAnyMap(self.dataDir + "/" + os.path.basename(mapfile))
			else:
				sys.exit("file " + mapfile + " not found")

		self.CheckUniqueFam()

        def GetPathList(self):
		return self.pathList;

        def GetPathName(self):
		return self.pathName;

	def ReadAnyMap(self, mapfile):
		try:
			file = open(mapfile, "r")
		except IOError:
			print "open file %s error" % mapfile

		for aline in file:
			if aline[0] == '#':
				continue
			subs = aline.strip().split()
			if len(subs) < 2:
				continue
			path, fam = subs[0], subs[1]
			if path in self.pathName:
				pathidx = self.pathName.index(path)
			else:
				pathidx = len(self.pathList)
				self.pathList.append(str(pathidx + 1))
				self.pathName.append(path)
				self.path2Fam.append([])
			if fam in self.famID:
				famidx = self.famID.index(fam)
			else:
				famidx = len(self.famList)
				self.famList.append(str(famidx + 1))
				self.famName.append(fam)
				self.famID.append(fam)
				self.fam2Path.append([])
			if famidx not in self.path2Fam[pathidx]:
				self.path2Fam[pathidx].append(famidx)
				self.fam2Path[famidx].append(pathidx)
				#same path & fun could be listed multiple times

		self.famTot = len(self.famList)
		self.pathTot = len(self.pathList)		
		print "total family", self.famTot, " pathway", self.pathTot
		#for idx in range(self.pathTot):
		#	print "path", self.pathName[idx], " include fam", len(self.path2Fam[idx])

	#read SEED figfam to subsytem mapping
	def ReadFigSubsytem(self, fig2ssfile):
		print "fig2ssfile=%s" % fig2ssfile
		try:
			file = open(fig2ssfile, "r")
		except IOError:
			print "open file %s error" % fig2ssfile
		for aline in file:
			aline = aline.strip()
			#note: in seed, a "family" (or "function") can have multiple FIG "sub"families
			#the families and the subsytems have names, but not codes
			#(subfam, fam, path) = aline.split('\t')
			cols = aline.split('\t')
			if len(cols) < 3: 
				continue  # skip the families that have NO subsystems assignment
			fam = cols[1]
			if fam in self.famName:
				famidx = self.famName.index(fam)	
			else:
				famidx = len(self.famName)
				self.famList.append("F" + str(famidx + 1))
				self.famID.append(cols[0])
				self.famName.append(fam)
				row = []
				self.fam2Path.append(row)
			for path in cols[2:]: 
				# note: a line could have multiple subsytems??
				if path[:4] == "CBSS":
					continue  # skip the clustering-based subsystems
				if path in self.pathName:
					pathidx = self.pathName.index(path)
				else:
					pathidx = len(self.pathName)
					self.pathList.append("S" + str(pathidx + 1))
					self.pathName.append(path)
					row = []
					self.path2Fam.append(row)

				if pathidx not in self.fam2Path[famidx]:
					self.fam2Path[famidx].append(pathidx)
				if famidx not in self.path2Fam[pathidx]:
					self.path2Fam[pathidx].append(famidx)
		file.close()
		self.pathTot = len(self.pathList)
		self.famTot = len(self.famList)
		print "total SEED subsystem=%d" % self.pathTot
		print "total SEED functions(families)=%d" % self.famTot

	#read KEGG pathways from ~/pathway/map-title.tab
	def ReadKEGGPath(self, pathfile):
		try:
			file = open(pathfile, "r")	
		except IOError:
			print "open file %s error" % pathfile
			sys.exit()

		for aline in file:
			aline = aline[:-1]
			row = aline.split("\t")
			self.pathList.append(row[0])
			self.pathName.append(row[1])
			row = []
			self.path2Fam.append(row)
		file.close()

		self.pathTot = len(self.pathList)
		print "total KEGG pathway=%d" % self.pathTot

	#read the KEGG families (KO) from ~/genes/ko
	def ReadKO(self, kofile = "", ifreadgene = False, ifgivenspe = ""):
		try:
			file = open(kofile, "r")	
		except IOError:
			print "open file %s error" % kofile
		print "read kofile=%s" % kofile
		for aline in file:
			aline = aline.strip()
			m = re.search(r'^ENTRY\s+(?P<ko>[\S]+)', aline)
			if m:
				self.famList.append(m.group('ko'))
				koidx = len(self.famList) - 1
				row = []
				self.fam2Path.append(row)
				for aline in file:
					name = aline[12:-1]
					self.famName.append(name)
					break
				#print "ko=%d %s name=%s" % (koidx, m.group('ko'), name)
				continue
			m = re.search(r'\[PATH:ko(?P<path>[^\]]+)', aline)
			if m:
				thispath = m.group('path')
				if thispath in self.pathList:
					pathidx = self.pathList.index(thispath)
					#print "thispath=%s %d idx=%d" % (thispath, pathidx, len(self.famList)) 
					self.fam2Path[koidx].append(pathidx)
					self.path2Fam[pathidx].append(koidx)
				continue

			#this information is only needed for KEGG pathway inference (blast-result to ko assignment)
			#reading this information is time-consuming
			if not ifreadgene:
				continue
       	 		m = re.search(r'^GENES', aline)
			if not m:
				continue
			#raw_input(" continue with reading gene...")
			lines = [aline]
			for aline in file:
				if aline[:3] == "///":
					break
				lines.append(aline)
			ifvalid = False
			for aline in lines:
                                #print 'now check line', aline
				if aline[15] == ':':
					org = aline[12:15].lower()
                                        #print "org", org, "ifgivenspe", ifgivenspe
					if ifgivenspe != "" and org != ifgivenspe:
						ifvalid = False 
					elif (ifgivenspe != "" and org == ifgivenspe) or (ifgivenspe == ""):
						ifvalid = True 
						if org in self.orgList:
							orgidx = self.orgList.index(org)
						else:
							orgidx = len(self.orgList)
							self.orgList.append(org)
							row = []
							self.orgGeneList.append(row)
							row = []
							self.orgGene2Fam.append(row)
				if not ifvalid:
					continue
				#else: the org is the same
                       		info = aline[17:]
                        	info2 = re.sub(r'\([^\)]+\)', "", info)
	                        items = info2.split()
        	                for aquery in items:
					#thisgene = org + ":" + aquery
					thisgene = aquery
					#print "thisgene = ", thisgene
					if thisgene in self.orgGeneList[orgidx]:
						geneidx = self.orgGeneList[orgidx].index(thisgene)
					else:
						geneidx = len(self.orgGeneList[orgidx])
						self.orgGeneList[orgidx].append(thisgene)
						row = []
						self.orgGene2Fam[orgidx].append(row)
					self.orgGene2Fam[orgidx][geneidx].append(koidx)

		file.close()
		self.famTot = len(self.famList)
		self.famID = self.famList
		self.orgTot = len(self.orgList)
		print "total KEGG fam=%d" % (self.famTot)
		if ifreadgene:
			print "total organisms involved =%d" % (self.orgTot)
			one2one = 0
			one2mul = 0
			for i in range(self.orgTot):
				#print " org %s has %d genes assigned to KO" % (self.orgList[i], len(self.orgGeneList[i]))
				for j in range(len(self.orgGeneList[i])):
					if len(self.orgGene2Fam[i][j]) > 1:
						one2mul += 1
					else:
						one2one += 1
			print "total %d genes matched to multiple KO; %d matched to single KO" % (one2mul, one2one)

	#calcualte the uniqueness of Fam (KEGG families or SEED families)
	def CheckUniqueFam(self):
		for i in range(self.famTot):
			if len(self.fam2Path[i]) > 10:
				print "fam=%d %s [%s] map-to-path=%d" % (i, self.famList[i], self.famName[i], len(self.fam2Path[i]))			
				for p in self.fam2Path[i]:
					print "   path-%d[%s %s]" % (p, self.pathList[p], self.pathName[p])
		ubifam = 0
		for i in range(self.famTot):
			if len(self.fam2Path[i]) > 1:
				ubifam += 1
		print "total families that mapped to more than one pathway = %d" % ubifam
					
		for p in range(self.pathTot):
			row = []
			totfam = len(self.path2Fam[p])
			uniquefam = 0 
			print "pathway-%d[%s; %s] fam=%d" % (p, self.pathList[p], self.pathName[p], totfam)
			for k in range(totfam):
				fam = self.path2Fam[p][k]
				if len(self.fam2Path[fam]) == 1:
					uniquefam = uniquefam + 1 
					row.append(fam)
					print "   unique fam-%d %s %s" % (fam, self.famList[fam], self.famName[fam])
			self.path2FamUni.append(row)
			print ">>>pathway-%d[%s; %s] fam=%d unique-fam=%d" % (p, self.pathList[p], self.pathName[p], totfam, uniquefam)
	def GetPath2FamUniMapped(self, apath, what):
		if apath in self.pathList:
			p = self.pathList.index(apath)	
			famlist = []
			for fam in self.path2FamUni[p]:
				if fam in self.famMapped:
					if what == "name":
						famlist.append(self.famList[fam])	
					else:
						famlist.append(fam)	
			return famlist
		else:
			return []

	def GetPath2FamMapped(self, apath, what):
		if apath in self.pathList:
			p = self.pathList.index(apath)	
			famlist = []
			#print "apath=", apath, "total fam", len(self.path2Fam[p])
			for fam in self.path2Fam[p]:
				#print "  >>>check fam", fam, "famMapped-total", len(self.famMapped)
				if fam in self.famMapped:
					if what == "name":
						famlist.append(self.famList[fam])
					else:
						famlist.append(fam)
					#print "   >>>Found in the famMapped"
			return famlist
		else:
			return []

	#return the index of families
	def GetPath2FamUni(self, apath):
		if apath in self.pathList:
			p = self.pathList.index(apath)	
			return self.path2FamUni[p]
		else:
			return []
	#return the index of families
	def GetPath2Fam(self, apath):
		if apath in self.pathList:
			p = self.pathList.index(apath)	
			return self.path2Fam[p]
		else:
			return []

	#note: this function only works when the givenspe is defined (see KEGG2html.py)
	def OrthMapBasedOnKO(self, spe=""):
		if spe == "":
			return
		if spe not in self.orgList:
			return
		self.famMapped = []
		s = self.orgList.index(spe) 
		print "spe", spe, "total gene", len(self.orgGeneList[s])
		for g in range(len(self.orgGeneList[s])):
			for fam in self.orgGene2Fam[s][g]:
				if fam not in self.famMapped:
					self.famMapped.append(fam)
		print "total family=", len(self.famList), "total mapped=", len(self.famMapped)
		#raw_input("type enter to continue")

	def OrthMap(self, famidxlist=[], famnamelist=[], famcount = []):
		if len(famidxlist) == 0 and len(famnamelist) == 0:
			return []
		if famcount:
			self.famCount = [0] * self.famTot
		#get only the orthologs that are mappped to a pathway
		if len(famidxlist) > 0:
			famlist = famidxlist
			famref = self.famList
		else:
			famlist = famnamelist
			famref = self.famName
		orthtot0 = len(famlist)
		orthtotfind = 0
		orthtotmap = 0 
		self.famMapped = []
		pathmap = [0] * self.pathTot
		for i in range(len(famlist)):
			orth = famlist[i]
			if orth in famref:
				idx = famref.index(orth)
				if famcount:
					self.famCount[idx] = famcount[i] 
				orthtotfind += 1
				if len(self.fam2Path[idx]) > 0:
					self.famMapped.append(idx) 
					orthtotmap += 1
					for path in self.fam2Path[idx]:
						pathmap[path] = 1
		print "original ortholog=%d found-in-the-fam-list=%d found-in-the-fam-mapped-to-pathway=%d" % (orthtot0, orthtotfind, orthtotmap)
		self.pathMapped = []
		for p in range(self.pathTot):
			if pathmap[p]:
				self.pathMapped.append(p)

	#assign orthologs to pathways
	#strategy 1: first assign unique ones -- then other orthologs 
	def Orth2PathUni(self, famidxlist=[], famnamelist=[], famcount=[]):
		self.OrthMap(famidxlist = famidxlist, famnamelist = famnamelist, famcount=famcount)
		orthtotmap = len(self.famMapped)
		#sort the orthologs based on their "uniqueness"
		fam2path = [] 
		for i in self.famMapped:
			if len(self.fam2Path[i]) == 0:
				continue
			fam2path.append((len(self.fam2Path[i]), i)) 
                #print "fam2path = ", len(fam2path)
                #raw_input("type enter to continue..")

		fam2pathsorted = sorted(fam2path, key=operator.itemgetter(0))
		pathsort = map(operator.itemgetter(0), fam2pathsorted)
		famsort = map(operator.itemgetter(1), fam2pathsorted)

		#for i in range(len(famsort)):
		#	print "fam0 %d fam %d [%s; %s] topath %d" % (i, famsort[i], self.famList[famsort[i]], self.famName[famsort[i]], pathsort[i])

		#pathfam & pathfam0: the number of fam families assigned to each pathway
		pathfam = [0] * self.pathTot
		pathfam0 = [0] * self.pathTot

		#pathfam0: the number of fam assigned to each pathway (using all multiple assignments)
		for fam in famsort:
			for p in self.fam2Path[fam]:
				pathfam0[p] = pathfam0[p] + 1

		#pathfam: the number of fam assigned to each pathway considering the "uniqueness" of fam to each pathway
		maxhit = pathsort[-1] 
		print "the maximum number of pathways a family is assigned to=%d" % maxhit
		hit = 1
		unassigned = len(famsort)
		beg = 0
		annpath = []
		famassign = [-1] * orthtotmap
		while (hit <= maxhit) and (unassigned > 0):
			#print "check hit=%d beg=%d" % (hit, beg)
			for k in range(beg, orthtotmap):
				if famassign[k] != -1:
					continue
				if pathsort[k] > hit:
					break 
				fam = famsort[k]
				maxsaturate = -1 
				maxsaturate_p = 0
				#print "check k=%d fam=%d %s %s" % (k, fam, self.famList[fam], self.famName[fam])
				for p in self.fam2Path[fam]:
					saturate = 1.0 * pathfam[p] / len(self.fam2Path[fam])
					if(saturate > maxsaturate):
						maxsaturate = saturate
						maxsaturate_p = p
				famassign[k] = maxsaturate_p
				pathfam[maxsaturate_p] = pathfam[maxsaturate_p] + 1
			beg = k

			unassigned = 0
			for k in range(orthtotmap):
				if famassign[k] == -1:
					unassigned = unassigned + 1
			#print "try hit=%d unassigned=%d (tot=%d)" % (hit, unassigned, orthtotmap)
			#raw_input()
			hit = hit + 1

		annpath0 = 0
		annpath = 0
		self.pathMappedOpt = []
		for k in range(self.pathTot):
			if pathfam0[k] != 0:
				annpath0 = annpath0 + 1
			if pathfam[k] != 0:
				annpath = annpath + 1
				self.pathMappedOpt.append(k)

		print "total pathway %d (%d) is found, compared to %d (%d)" % (annpath, len(self.pathMappedOpt), annpath0, len(self.pathMapped))
		print "%-50s %s\t%s" % ("#pathway", "fam-assigned(all)", "fam-assigned(unique)[weight]")
		for p in range(self.pathTot):
			if pathfam0[p] == 0 and pathfam[p] == 0:
				continue
			tmp = self.pathList[p] + "[" + self.pathName[p] + "]"
			weight = .0
			for k in range(orthtotmap):
				ks = famsort[k]
				if famassign[k] != p:
					continue
				weight = weight + 1.0 / len(self.fam2Path[ks])
			print "%-50s %d\t%d[%.1f]" % (tmp, pathfam0[p], pathfam[p], weight)
		return self.pathMappedOpt

	#Parsinomy approach to pathway inference
	def Orth2PathMin(self, famidxlist=[], famnamelist=[], famcount=[], mpsfile="test.mps", glpsol=""):
		# write mps file (the input for glpsol, the integer programming package)
		print "now write mps file.."
		self.WriteMps(famidxlist=famidxlist, famnamelist=famnamelist, famcount=famcount, mpsfile=mpsfile)

		# run glpsol
		global glpsol0
		if glpsol == "":
			glpsol = glpsol0
		lpout = mpsfile + ".LPout"
		command = glpsol + " " + mpsfile + " -o " + lpout
		print "now run command = %s" % command
		os.system(command)

		# check the result
		self.GetLPOut(lpout)

		return self.pathMappedOpt

	#output mps file for integer programming (most parsinomy pathway inference)
	def WriteMps(self, famidxlist=[], famnamelist=[], famcount=[], mpsfile="test.mps"):
		try:
			file = open(mpsfile, "w")
		except IOError:
			print "open file %s error" % mpsfile
			sys.exit()
	        str = "%-14s%s\n" % ("NAME", "PATH")
       	 	file.write(str)
	
		self.OrthMap(famidxlist=famidxlist, famnamelist=famnamelist, famcount=famcount)
		orthtotmap = len(self.famMapped)

        	#write ROWS
		file.write("ROWS\n")
		file.write(" N  NUM\n")
		for orth in self.famMapped:
                	str = " G  F%s\n" % self.famList[orth];
			#note 1: use a different idx of family for mps file
                	#note 2: when use E, there is no feasible solution
                	#G>=1, mean each family has to be assigned to at least one pathway 
			file.write(str)

        	#write COLUMNS
		#the same column (pathway) needs to be organized in the same block 
        	file.write("COLUMNS\n")
		pathvalid = [0] * self.pathTot
		for p in self.pathMapped:
			#note: use a different pathway idx in mps
			pathname = "P%s" % self.pathList[p]
		        str = "    %-10s%-10s%10d\n" %(pathname, "NUM", 1)
			file.write(str)
			for orth in self.famMapped:
       	                 	famname = "F%s" % self.famList[orth];
				for path in self.fam2Path[orth]:
					if path == p:
		                        	str = "    %-10s%-10s%10d\n" %(pathname, famname, 1)
						file.write(str)	
						pathvalid[p] = 1

        	#write RHS
		file.write("RHS\n")
		for orth in self.famMapped:
                	famname = "F%s" % self.famList[orth]
                	str = "    %-10s%-10s%10.1f\n" % ("RHS1", famname, 1.0);
			file.write(str)

        	#write bounds
		file.write("BOUNDS\n")
		for p in range(self.pathTot):
			if pathvalid[p]:
				path = self.pathList[p]
	                	pathname = "P%s" % path;
       	         		str = " BV %-10s%-10s\n" %("BND1", pathname);
                		#all variants are binary (1 keep the pathway; 0 pathway not necessary)
				file.write(str)

        	file.write("ENDATA\n")
		file.close()		

        	print "End of PrintMPS"
		
	def GetLPOut(self, lpoutfile="test.mps.LPout"):
		try:
			file = open(lpoutfile, "r")
		except IOError:
			print "open file %s error" % lpoutfile 
			sys.exit(1)

		keeppath = []
		for aline in file:
			aline = aline.strip()
			cols = aline.split()
			if len(cols) < 2:
				continue
			if cols[0] == "Columns:":
				Columns = int(cols[1])
			elif cols[0] == "Objective:":
				MINimum = int(cols[3])
			elif cols[0] == "No." and cols[1] == "Column":	
				for aline2 in file:
					if aline2[0] == '-': 
						continue
					aline2 = aline2.strip()
					cols2 = aline2.split()
					if len(cols2) < 1:
						break
					if cols2[3] == '1':
						keeppath.append(cols2[1][1:])
						#check with WriteMps: the pathway idx used in mps is to add "P" before the pathList
		file.close()
		if len(keeppath) != MINimum:
			print "reading %s error: minimum %d read %d" % (lpoutfile, MINimum, len(keeppath)) 
			sys.exit()

		self.pathMappedOpt = []
		for path in keeppath:
			if path not in self.pathList:
				print "Error: unknown pathList %s" % path
				sys.exit()
			pathidx = self.pathList.index(path)
			self.pathMappedOpt.append(pathidx)	

		print "total pathways mappd: before inference %d, after inference %d" % (len(self.pathMapped), len(self.pathMappedOpt))

	#add the pathways with many functions annotated, even they were considered as redundant ones!
	def PopulatePath(self, pathmapped = [], par = 0.7):
		famvalid = [0] * self.famTot
		for fam in self.famMapped:
			famvalid[fam] = 1

		addpath = 0
		for p in range(self.pathTot):
			if len(self.path2Fam[p]) == 0:
				continue
			if not p in pathmapped:
				add = 0
				for f in self.path2Fam[p]:
					if famvalid[f] == 1:
						add += 1
				#pathways with most functions annotated should be added back -- even it is a redundant one
				print "pathway", p, self.pathList[p], self.pathName[p], "path2fam", len(self.path2Fam[p]), " real-family", add
				if add >= len(self.path2Fam[p]) * par:
					pathmapped.append(p)
					addpath += 1
					print "this pathway is added back"
				#else:
				#	print "this pathway does not have enough functions"
				#raw_input("type enter to continue")

		print "added pathway =", addpath

	#remove the pathways with too few functions annotated (e.g., 2, use par), 
	#even when their associated families are annotated(but NOT the unique ones)
        #not unique families assigned to this pathway? ubiquitous families have to be assigned to at least one of the pathways, right
	def RemoveSparsePath(self, pathmapped = [], par = 2):			
		famvalid = [0] * self.famTot
		for fam in self.famMapped:
			famvalid[fam] = 1

		delpath = 0
		for p in range(self.pathTot):
			if len(self.path2Fam[p]) == 0:
				continue
			if p in pathmapped:
				add = 0
				uni = 0
				for f in self.path2Fam[p]:
					if famvalid[f] == 1:
						add += 1
						if len(self.fam2Path[f]) == 1: 
							uni += 1
				#pathways with few functions annotated are removed
				if uni == 0 and add <= par:
					pathmapped.remove(p)
					delpath += 1
					print "pathway", p, self.pathList[p], self.pathName[p], "path2fam", len(self.path2Fam[p]), " real-family", add, " is removed from the list!!"
					#raw_input("type enter to continue")

		print "deleted pathway =", delpath

	def DiffPathMap(self, maps, tags):
		maps.insert(0, self.pathMapped)
		tags.insert(0, "Ori")
		mapnum = len(maps)
		pathvalid = intmatrix(self.pathTot, mapnum)
		for m in range(mapnum):
			for p in maps[m]:
				pathvalid[p][m] = 1
		print "#Summary for the pathway inference"
		#print description line
		str = "%-5s %-10s %-70s %-5s %-5s" % ("ID", "List", "Name", "Fam", "Fam-found")
		for atag in tags:
			str += " %-3s" % atag
		print str + " Same/Diff"
		#print each pathway
		totsame = 0
		totdiff = 0
		famvalid = [0] * self.famTot
		for fam in self.famMapped:
			famvalid[fam] = 1

		for p in range(self.pathTot):
			add = sum(pathvalid[p])
			if add == 0:
				continue
			add = sum(pathvalid[p][1:])
			if add == mapnum - 1 or add == 0:
				label = "Same" 
				if add != 0:
					totsame += 1
			else:
				label = "Diff"
				totdiff += 1
			add = 0
			for f in self.path2Fam[p]:
				if famvalid[f] == 1:
					add += 1
			str = "%-5d %-10s %-70s %-5d %-5d" % (p + 1, self.pathList[p], self.pathName[p], len(self.path2Fam[p]), add)
			for m in range(len(maps)):
				str += " %-3d" % pathvalid[p][m] 
			print str + " " + label
		#print total number line
		str = "%-5s %-10s %-70s %-5s %-5s" % ("#total", "", "", "", "")
		for m in range(mapnum):
			str += " %-3d" % len(maps[m]) 
		print str
		#print functional diveristy line
		#str = "#functional-diversity [max: log(%d)=%.3f]" % (self.pathTot, math.log(self.pathTot * 1.0))
		#str = "%-99s" % str
		#for m in range(mapnum):
		##	str += " %.3f" % math.log(len(maps[m]) * 1.0)
		#print str
		print "#total match=%d diff=%d" % (totsame, totdiff)

	def WriteReport(self, minpath, reportfile, detailfile):
		na = True
		if self.whichDB == "KEGG":
			keggmap = [] 
			keggdir = "/dataomics/kegg/kegg-curr"
			mapfile = keggdir + "/pathway/" + self.speID.lower() + "/map.list"
			if os.path.exists(mapfile):
				file = open(mapfile, "r")
				for aline in file:
					m = re.match('^[^\d]+(?P<id>\d+)', aline)
					if m:
						id = m.group('id')
						idx = self.pathList.index(id)
						keggmap.append(idx)
				na = False
			tags = ["kegg", "naive", "minpath"]
			maps = [keggmap, self.pathMapped, minpath]
		else:
			seedmap = []	
			if self.whichDB == 'SEED':
				tags = ["seed", "naive", "minpath"]
			else:
				tags = ["any", "naive", "minpath"]
			maps = [seedmap, self.pathMapped, minpath]
		mapnum = len(maps)
		pathvalid = intmatrix(self.pathTot, mapnum)
		for m in range(mapnum):
			for p in maps[m]:
				pathvalid[p][m] = 1
		famvalid = [0] * self.famTot
		for fam in self.famMapped:
			famvalid[fam] = 1

		file = open(reportfile, "w")
		if detailfile:
			detail = open(detailfile, "w")
		for p in range(self.pathTot):
			add = sum(pathvalid[p])
			if add == 0:
				continue
			add = 0
			for f in self.path2Fam[p]:
				if famvalid[f] == 1:
					add += 1
			if na:
				tmp = "n/a"
			else:
				tmp = pathvalid[p][0]
			print >> file, "path", self.pathList[p], tags[0], tmp, " naive", pathvalid[p][1], " minpath", pathvalid[p][2], " fam0 ", len(self.path2Fam[p]), " fam-found ", add, " name ", self.pathName[p]

			if not (detailfile and pathvalid[p][2]):
				continue
			#print details
			print >> detail, "path", self.pathList[p], "fam0", len(self.path2Fam[p]), "fam-found", add, "#", self.pathName[p]
			for f in self.path2Fam[p]:
				if famvalid[f] == 1 and self.famCount:
					print >> detail, "  ", self.famID[f], "hits", self.famCount[f], "#", self.famName[f]
				elif famvalid[f] == 1:
					print >> detail, "  ", self.famID[f], "#", self.famName[f]
					
		print
		print "Results are saved in file:", reportfile 
		file.close()
		if detailfile:
			print "Details are saved in file:", detailfile
			detail.close()

# this function reads in KO/fig assignment, then map KO/fig families to the pathways
# last update by Yuzhen Ye on July 3, 2009
def Orth2Path(infile = "demo.ko", whichdb = "KEGG", mpsfile = "test.mps", reportfile = "test.minpath", detailfile = "", mapfile=""):
	try:
		file = open(infile, "r")
	except IOError:
		sys.exit( "open file error " + infile)

	orthlist, orthcount = [], []
	add = 0
	for aline in file:
		tmp = aline.strip().split("\t") #Jan 11, 2018
		if len(tmp) < 2:
			tmp = aline.split()
		if len(tmp) >= 2:
			add = add + 1
			if tmp[1] not in orthlist:
				orthlist.append(tmp[1])
				orthcount.append(1)
				#print "ko-%d=%s" % (len(orthlist), tmp[1])
			else:
				idx = orthlist.index(tmp[1])
				orthcount[idx] += 1
	file.close()

	print "total input orth=%d  unique=%d" % (add, len(orthlist))

	test = MinPath(whichdb = whichdb, mapfile = mapfile) #default pathwaydb: KEGG

	if whichdb == "KEGG":
		map = test.Orth2PathMin(famidxlist=orthlist, famnamelist=[], famcount=orthcount, mpsfile=mpsfile)
	elif whichdb == "SEED":
		map = test.Orth2PathMin(famidxlist=[], famnamelist=orthlist, famcount=orthcount, mpsfile=mpsfile)
	else:
		map = test.Orth2PathMin(famidxlist=[], famnamelist=orthlist, famcount=orthcount, mpsfile=mpsfile)
	#KEGG by ids, and fig by names
		
	map_add = map[:]
	par = 0.5
	test.PopulatePath(pathmapped = map_add, par = par)

	test.WriteReport(map_add, reportfile, detailfile)

	os.system("rm test.mps*")

if __name__ == '__main__':
	kofile, figfile, anyfile, mapfile, mpsfile, reportfile, detailfile = "", "", "", "", "test.mps", "test.minpath", ""
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-ko":
			kofile = sys.argv[i + 1]
		elif sys.argv[i] == "-fig":
			figfile = sys.argv[i + 1]
		elif sys.argv[i] == "-any":
			anyfile = sys.argv[i + 1]
		elif sys.argv[i] == "-map":
			mapfile = sys.argv[i + 1]
		elif sys.argv[i] == "-report":
			reportfile = sys.argv[i + 1]
		elif sys.argv[i] == "-details":
			detailfile = sys.argv[i + 1]
		elif sys.argv[i] == "-mps":
			mpsfile = sys.argv[i + 1]
	if kofile:
		Orth2Path(infile = kofile, mpsfile = mpsfile, reportfile = reportfile, detailfile = detailfile)
	elif figfile:
		Orth2Path(infile = figfile, mpsfile = mpsfile, reportfile = reportfile, detailfile = detailfile, whichdb = "SEED")
	elif anyfile and mapfile:
		Orth2Path(infile = anyfile, mpsfile = mpsfile, reportfile = reportfile, detailfile = detailfile, whichdb = "ANY", mapfile=mapfile)
	else:
		print "Usage: python MinPath.py <-ko filename>/<-fig filename>/<-any annfile> [-map mapfile] [-report filename] [-details detailed-output]"
		print "Note: your input file can contain functional annotations in either of the following"
		print "   -ko file: annotation in KEGG KO families"
		print "   -fig file: annotation in SEED fig families"
		print "   -any file: annotation in any families, then you must specify -map, the pathway-function mapping file"
		print "Example 1: python MinPath.py -ko demo.ko -report demo.ko.minpath"
		print "Example 2: python MinPath.py -ko demo.ko -report demo.ko.minpath -details demo.ko.minpath.details"
		print "Example 3: python MinPath.py -fig demo.fig -report demo.fig.minpath"
		print "Example 4: python MinPath.py -fig demo.fig -report demo.fig.minpath -details demo.fig.minpath.details"
		print "Example 5: python MinPath.py -any demo.ec -map ec2path -report demo.ec.minpath -details demo.ec.minpath.details"
		sys.exit(1)

