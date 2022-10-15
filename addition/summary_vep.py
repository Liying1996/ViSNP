# Input: results generated from VEP
# Output: Summary file for following R code
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", help = "VEP annotation", required = True)
parser.add_argument("-o", help = "Output file name", required = True)

args = parser.parse_args()

input_file = args.i
output_file = args.o

titles = []
new = open(output_file, "w")

flag = False
with open(input_file) as f:
	for line in f:
		line = line.strip()
		if line[:2] == "##":
			title = line.split(':')[0].replace("##", "").strip()

			if title == "Uploaded_variation":
				flag = True
			
			if flag:
				if title != "Extra column keys":	
					titles.append(title)

		elif line[:2] == "#U":
			new.write("\t".join(titles) + "\n")

		else:
			front_line = line.split("\t")[:-1]
			write_line = "\t".join(front_line)

			extras_line = line.split("\t")[-1].split(";")
			
			info = dict()
			for each in extras_line:
				name = each.split("=")[0]
				content = each.split("=")[-1]
				info[name] = content

			for t in titles[13:]:
				if t in info:
					write_line += "\t" + info[t]
				else:
					write_line += "\t-"

			new.write(write_line + "\n")
new.close()
