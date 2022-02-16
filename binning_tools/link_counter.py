import re
import sys
args = sys.argv

def make_scaffold_list(samfile,minimumlength):
	sclist = {}
	with open(samfile,"r")as data:
		for line in data:
			line_trim = line.strip("\n")
			table = line_trim.split("\t")
			if table[0] == "@SQ":
				if int(table[2][3:]) >= minimumlength:
					sclist.setdefault(int(re.search('[0-9]{1,8}',table[1]).group()),int(table[2][3:]))
			if line[0] != "@":
				break
	return sclist

def make_contact_number(samfile,sclist):
	contact_number = {}
	with open(samfile,"r")as data:
		switch = 0
		read = []
		for line in data:
			if line[0] != "@":
				line_trim = line.strip("\n")
				table = line_trim.split("\t")
				if switch == 0:
					read.append(table)
					switch = 1
				if switch == 1:
					if table[0] == read[0][0]:
						read.append(table)
						continue
					if table[0] != read[0][0]:
						if len(read) != 2:
							read = [table]
							continue
						if len(read) == 2:
							forward = int(re.search('[0-9]{1,8}',read[0][2]).group())
							reverse = int(re.search('[0-9]{1,8}',read[1][2]).group())
							if forward in sclist and reverse in sclist:
								if forward != reverse:
									if int(read[0][4]) == 60 and int(read[1][4]) == 60:
										key = ""
										if forward >= reverse:
											key = str(reverse) + "-" + str(forward)
										if forward < reverse:
											key = str(forward) + "-" + str(reverse)
										if key in contact_number:
											contact_number[key] += 1
										contact_number.setdefault(key,1)
						read = [table]
		if len(read) == 2:
			forward = int(re.search('[0-9]{1,8}',read[0][2]).group())
			reverse = int(re.search('[0-9]{1,8}',read[1][2]).group())
			if forward in sclist and reverse in sclist:
				if forward != reverse:
					if int(read[0][4]) == 60 and int(read[1][4]) == 60:
						key = ""
						if forward >= reverse:
							key = str(reverse) + "-" + str(forward)
						if forward < reverse:
							key = str(forward) + "-" + str(reverse)
						if key in contact_number:
							contact_number[key] += 1
						contact_number.setdefault(key,1)
	return contact_number

def make_sorted_dict(contact_number):
	contact_number_sorted = sorted(contact_number.items(),key=lambda x:-x[1])
	contact_number_dict = {}
	for i in range(len(contact_number_sorted)):
		contact_number_dict.setdefault(contact_number_sorted[i][0],contact_number_sorted[i][1])
	edgelist = {}
	label_edge = {}
	edge_number = 0
	for k,v in contact_number_dict.items():
		label = k.split("-")
		source = int(label[0])
		target = int(label[1])
		if source in label_edge:
			label_edge[source].append(edge_number)
		label_edge.setdefault(source,[edge_number])
		if target in label_edge:
			label_edge[target].append(edge_number)
		label_edge.setdefault(target,[edge_number])
		edgelist.setdefault(edge_number,[source,target,v])
		edge_number += 1
	edge_selected = {}
	for k,v in label_edge.items():
		if sclist[k] >= 300000:
			for i in range(len(v)):
				if i == 0:
					if edgelist[v[i]][2] >= 50:
						if v[i] in edge_selected:
							edge_selected[v[i]] += 1
						edge_selected.setdefault(v[i],1)
				if i == 1:
					if edgelist[v[i]][2] >= 50 and (edgelist[v[i]][2]/edgelist[v[0]][2]) >= 0.2:
						if v[i] in edge_selected:
							edge_selected[v[i]] += 1
						edge_selected.setdefault(v[i],1)
				if i > 1:
					if edgelist[v[i]][2] >= 50 and (edgelist[v[i]][2]/edgelist[v[0]][2]) >= 0.25:
						if v[i] in edge_selected:
							edge_selected[v[i]] += 1
						edge_selected.setdefault(v[i],1)
		if sclist[k] < 300000:
			for i in range(len(v)):
				if i == 0:
					if edgelist[v[i]][2] >= 5:
						if v[i] in edge_selected:
							edge_selected[v[i]] += 1
						edge_selected.setdefault(v[i],1)
				if i == 1:
					if edgelist[v[i]][2] >= 10:
						if v[i] in edge_selected:
							edge_selected[v[i]] += 1
						edge_selected.setdefault(v[i],1)
				if i > 1:
					if edgelist[v[i]][2] >= 15:
						if v[i] in edge_selected:
							edge_selected[v[i]] += 1
						edge_selected.setdefault(v[i],1)
	contact_number_selected = {}
	for k,v in edge_selected.items():
		if v == 2:
			key = str(edgelist[k][0]) + "-" + str(edgelist[k][1])
			contact_number_selected.setdefault(key,edgelist[k][2])
	return contact_number_selected

def making_scdict(sc_fa):
	scdict = {}
	with open(sc_fa,"r")as data:
		temporary_seq = ""
		scnum = 0
		for line in data:
			line_trim = line.strip("\n")
			if line_trim[0] != ">":
				temporary_seq += line_trim
			if line_trim[0] == ">":
				nodenum = int(re.search("[0-9]{1,8}",line_trim).group())
				if temporary_seq == "":
					scnum = nodenum
				if temporary_seq != "":
					scdict.setdefault(scnum,temporary_seq)
					scnum = nodenum
				temporary_seq = ""
		if temporary_seq != "":
			scdict.setdefault(scnum,temporary_seq)
	return scdict

def search_cutsite(scdict,cuttingsite):
	sc_site = {}
	count = 0
	for k,v in scdict.items():
		for i in range(len(v)-len(cuttingsite)+1):
			if cuttingsite == v[i:i+len(cuttingsite)]:
				count += 1
		if count == 0:
			count = 1
		sc_site.setdefault(k,count)
	return sc_site

def correct(contact_number,sc_site):
	sc_corrected = {}
	for k,v in contact_number.items():
		table = k.split("-")
		for i in range(len(table)):
			table[i] = int(table[i])
		correction = sc_site[table[0]] * sc_site[table[1]]
		v_corrected = 0.0
		v_corrected = v/correction
		sc_corrected.setdefault(k,v_corrected)
	return sc_corrected

def make_contactfile(sc_corrected,contactfile):
	with open(contactfile,"w")as f:
		for k,v in sc_corrected.items():
			line = str(k) + "\t" + str(v) + "\n"
			f.write(line)

sclist = make_scaffold_list(args[6],int(args[5]))
contact_number = make_contact_number(args[1],sclist)
contact_number_sorted = make_sorted_dict(contact_number)
scdict = making_scdict(args[2])
sc_site = search_cutsite(scdict,args[3])
sc_corrected = correct(contact_number_sorted,sc_site)
make_contactfile(sc_corrected,args[4])
"""
args[1] -> samfile name
args[2] -> fastafile name
args[3] -> cutting_site
args[4] -> output file name
args[5] -> mimimum length
args[6] -> header.txt
"""
