import matplotlib.pyplot as plt
import json
import numpy as np
import pandas as pd
MAX_FEATURE = 576
PARENT_FOLDER = '/Users/dcongtinh/gene-abundance/experiment/results/fc_model_with_feature_selected/qtf_pc576_10_fillseqf_nb10_auy_gray/'

def read(PATH):
	f=open(PATH, 'r')
	freq = {}
	for line in f.readlines():
		obj = json.loads(line)
		for _ in obj:
			if _ in freq:
				freq[_]+=1
			else:
				freq[_] = 1
	return freq

def export2csv(DATA, EXPORT_PATH, DATASET_NAME):
	freq = {'freq': DATA}
	df = pd.DataFrame(freq, columns= ['freq'])
	export_csv = df.to_csv (EXPORT_PATH + DATASET_NAME, index=False, header=True)
	print ("Data was exported at %s" % (EXPORT_PATH))

def main():
	for dataset in ['cirgene', 'colgene', 'ibdgene', 'obegene', 't2dgene', 'wt2dgene']:
		freq = read(PARENT_FOLDER + 'feature_selected_' + dataset + '.txt')
		sorted = []
		for key in freq: 
			sorted.append((freq[key], key))
		sorted.sort(reverse=True)
		sorted = sorted[:MAX_FEATURE]
		__freq__ = []
		for value, key in sorted:
			__freq__.append(value)
		export2csv(__freq__, PARENT_FOLDER, 'feature_selected_' + dataset + '.csv')
	# top = sorted[:MAX_FEATURE]

	# plt.bar([i for i in range(len(top))], [value for value, key in top], tick_label=[key for value, key in top])
	# plt.savefig(INPUT_FILE+'.png')
	# plt.show()
	pass
if __name__=='__main__':
	main()	