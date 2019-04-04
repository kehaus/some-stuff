

import os


def redact(line):
	return line.replace('\t','    ').replace('\{','{').replace('\}','}').replace('\\','').replace('{textbackslash}','\\')
    #return line.replace('\t','    ').replace('{\\textbackslash}','\\').replace('\{','{').replace('\}','}').replace('\\','')

list_ = open('PennState.bib').read().split('\n')



list_2 = [redact(line) for line in list_]

print('test')
print(list_2)
print('test')

with open('zzz3.bib', 'a') as f:
	for line in list_2:
		f.write(line)
		f.write('\n')

print('test')
