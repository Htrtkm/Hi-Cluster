import sys
import re

def convert_contactfile(contactfile):
	contact_number = {}
	with open(contactfile,"r")as data:
		for line in data:
			line_trim = line.strip("\n")
			table = line_trim.split()
			contact_number.setdefault(table[0],table[1])
	return contact_number

def merge_contact_number(contactnumber1,contactnumber2):
	merged_contact_number = {}
	for k,v in contactnumber1.items():
		merged_contact_number.setdefault(k,float(v))
	for k,v in contactnumber2.items():
		if k in merged_contact_number:
			merged_contact_number[k] += float(v)
		merged_contact_number.setdefault(k,float(v))
	return merged_contact_number

def making_contig_sorted(merged_contact_number):
	contignum = {}
	for k,v in merged_contact_number.items():
		table = k.split("-")
		X = int(table[0])
		contignum.setdefault(X,[])
		contignum[X].append([table[1],v])
	contact_number_sorted = sorted(contignum.items(),key = lambda x:x[0])
	contignum = {}
	return contact_number_sorted

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

def making_gmlfile(inputfile_sorted,outputfile,sclist):
	with open(outputfile,"w")as f:
		nodeinfo = ""
		edgeinfo = ""
		f.write("Creator \"Takumi Hattori\"" + "\n" + "graph" + "\n" + "[" + "\n")
		node_list = {}
		for i in range(len(inputfile_sorted)):
			A = inputfile_sorted[i][0]
			B = inputfile_sorted[i][1]
			node_list.setdefault(int(A),0)
			for j in range(len(B)):
				node_list.setdefault(int(B[j][0]),0)
		node_list_sorted = sorted(node_list.items(),key = lambda x:x[0])
		id_num = 0
		for i in range(len(node_list_sorted)):
			label = "label \"" + str(node_list_sorted[i][0]) + "\""
			nodeinfo += "\t" + "node" + "\n" + "\t" + "[" + "\n" + "\t" + "\t" + "id " + str(i) + "\n" + "\t" + "\t" + label + "\n" + "\t" + "]" + "\n"
			node_list[node_list_sorted[i][0]] = i
			id_num = i
		for k,v in sclist.items():
			if k not in node_list:
				id_num += 1
				nodeinfo += "\t" + "node" + "\n" + "\t" + "[" + "\n" + "\t" + "\t" + "id " + str(id_num) + "\n" + "\t" + "\t" + "label \"" + str(k) + "\"" + "\n" + "\t" + "]" + "\n"
		for i in range(len(inputfile_sorted)):
			A = inputfile_sorted[i][0]
			B = inputfile_sorted[i][1]
			for j in range(len(B)):
				value = B[j][1]
				#value = math.log10(value)
				#normalized_value = (value / (cut_site_list[int(A)] * cut_site_list[int(B[j][0])])) * (10**10)
				#value not normalized
				source_target = "\t" + "\t" + "source " + str(node_list[int(A)]) + "\n" +"\t" +  "\t" + "target " + str(node_list[int(B[j][0])]) + "\n" + "\t" + "\t" + "weight " + str(value)
				edgeinfo += "\t" + "edge" + "\n" + "\t" + "[" + "\n" + source_target + "\n" + "\t" + "]" + "\n"
		f.write(nodeinfo)
		f.write(edgeinfo)
		f.write("]")

args = sys.argv

sclist = make_scaffold_list(args[4],int(args[5]))
contact_num1 = convert_contactfile(args[1])
contact_num2 = convert_contactfile(args[2])
merged_contact_number = merge_contact_number(contact_num1,contact_num2)
contact_number_sorted = making_contig_sorted(merged_contact_number)
making_gmlfile(contact_number_sorted,args[3],sclist)
"""
args[1] => contact file
args[2] => contact file
args[3] => output file
args[4] => header.txt
args[5] => minimum length
"""
