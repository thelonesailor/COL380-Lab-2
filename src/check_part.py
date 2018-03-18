import sys

if(len(sys.argv) != 2):
	print("Usage: python check_part.py <output>")
	exit()

part_map = {}

fin = open(sys.argv[1], 'r')

line = fin.readline()

lsp = line.split()

tot_num = len(lsp)
for w in lsp:
	num = int(w)
	if(num in part_map):
		part_map[num] += 1
	else:
		part_map[num] = 1

print("Total vertices = {}".format(tot_num))

l = []
print("All class counts:")
for x in part_map:
	print("Class {} has count {}".format(x, part_map[x]))
	l.append(part_map[x])

x_max = max(l); x_min = min(l)

thresh = abs(float(x_min) - float(x_max))*1.0/float(tot_num)

if(thresh <= 0.05):
	print("PASS")
else:
	print("FAIL")