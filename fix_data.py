print("started")
## read in file, split elements of each row into two lists
f = open("african_readlines2.txt", "r+")
new_str = ""
str_list = []
string = []
build = f.readlines()
for x in range(len(build)-1):
    for y in range(len(build[x])):
        if build[x][y] != "	" and build[x][y] != "\t":
            new_str += build[x][y]
        if build[x][y] == "	" or build[x][y] == "\t":
            new_str += "\t"
            string.append(new_str)
            new_str = ""
    str_list.append(string)
    string = []

mat_list = []
pat_list = []
id = ""
for t in range(int(len(str_list)/2)-1):
    if len(str_list[t]) != 0:
        id = str(str_list[t][0]) + str(str_list[t][1])
        for y in range(3, len(str_list[t])):
            if y % 2 != 0:
                mat_list.append(str_list[t][y])
            else:
                pat_list.append(str_list[t][y])
        g = open("african_fixed5_2.txt", "a+")
        g.write(str(id))
        for z in range(len(mat_list)-1):
            g.write(str(mat_list[z]))
        g.write(str(id))
        for a in range(len(pat_list)-1):
            g.write(str(pat_list[a]))
        mat_list = []
        pat_list = []
        print("element " + str(t) + " successfully copied into two individuals")

g = open("african_fixed5_2.txt", "r+")
h = open("african_fixed6_2.txt", "w+")
lines = g.readlines()
for x in range(len(lines)-1):
    if len(lines[x]) != 0 and lines[x] != "\n" and lines[x] != None:
        if lines[x][0] == '-':
            if lines[x][1] == '9':
                lines[x] = lines[x][2:len(lines[x])-1]
        else:
            h.write(lines[x])

print("done!")