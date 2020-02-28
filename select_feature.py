import matplotlib.pyplot as plt
import json
import numpy as np
MAX_FEATURE = 576
INPUT_FILE = 'sample.txt'

def read():
	f=open(INPUT_FILE, 'r')
	freq = {}
	for line in f.readlines():
		obj = json.loads(line)
		for _ in obj:
			if _ in freq:
				freq[_]+=1
			else:
				freq[_] = 1
	return freq
def main():
	freq = read()
	sorted = []
	for key in freq: 
		sorted.append((freq[key], key))
	sorted.sort(reverse=True)
	top = sorted[:MAX_FEATURE]
	print(top)
	plt.bar([i for i in range(len(top))], [value for value, key in top], tick_label=[key for value, key in top])
	plt.savefig(INPUT_FILE+'.png')
	plt.show()
	pass
if __name__=='__main__':
	main()	