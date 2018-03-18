import sys

if(len(sys.argv) != 2):
	print("Usage: python check_part.py <output>")
	exit()

part_map = {}

fin = open(sys.agrv[1], 'r')

line = f.readline()

lsp = line.split()

tot_num = len(slp)
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

thresh = abs(double(x_min) - double(x_max))*1.0/double(tot_num)

if(thresh <= 0.05):
	print("PASS")
else:
	print("FAIL")