'''
Use Deepmg
Created: 21/02/2020
'''

import os
import time

dataset = ['cirgene', 'colgene', 'ibdgene', 'obegene', 't2dgene', 'wt2dgene'] # Gene family datasets

# Set the PATH:
ROOT = '/Users/dcongtinh'
path_data = ROOT + '/Downloads/data/'                                   # For my macbook
machine_path = ROOT + '/gene-abundance'                                 # For my macbook
file_run = ROOT +'/gene-abundance/subdeepmg/__main__.py'                # Use package deepmg to modify
file_read_results = ROOT +'/gene-abundance/subdeepmg/read_results.py'   # Use package deepmg to modify
algo_redu = 'pwbf' # Perceptron Weight Based Filter
new_dim = 576   # Image 32x32 pixels
run_time = 10
num_bin = 10
type_bin = 'eqf'
scale_mode = 'qtf'
models = ['fc_model', 'model_cnn']
colormaps = ['gray', 'rainbow']
type_emb = 'fills'

cnt = 0
def executeCommand(__command__):
    global cnt
    cnt += 1
    start_time = time.time()
    os.system(__command__)
    end_time = time.time()
    print("Run-time: %f seconds"%(end_time-start_time))
    os.system('say Done %d' % cnt)
    # print(__command__ + '\n')

start_time = time.time()
time_start = time.ctime()
for model_v in models:
    for colormap in colormaps:
        for dataset_name in dataset:
            res_folder = '%s/experiment/results/%s/%s_%s%d_10_%s%s_nb%d_auy_%s/' % (machine_path, model_v, scale_mode, algo_redu, new_dim, type_emb, type_bin, num_bin, colormap)   # Result folder
            img_folder = '%s/experiment/images/%s/%s_%s%d_10_%s%s_nb%d_auy_%s/' % (machine_path, model_v, scale_mode, algo_redu, new_dim, type_emb, type_bin, num_bin, colormap)    # Image folder

            command = 'python3 %s --save_para y --original_data_folder %s --data_name %s --run_time %d --algo_redu %s --new_dim %d --model %s --auto_v y --parent_folder_results %s --parent_folder_img %s --scale_mode %s --type_bin %s --num_bin %d --type_emb %s -z 255 --colormap %s' % (file_run, path_data, dataset_name, run_time, algo_redu, new_dim, model_v, res_folder, img_folder, scale_mode, type_bin, num_bin, type_emb, colormap)
            if colormap == 'gray':      # Gray image
                command += ' --channel 1'
            else:                       # Color image
                command += ' --channel 3 --preprocess_img vgg16'
            executeCommand(command)
            os.system('python3 %s -i %s -o %s' % (file_read_results, res_folder, res_folder))
            exit()

print("Training dataset use [%s] at %s\n.\n.\n." % (algo_redu.upper(), time_start))
print("Totsl run-time: %f seconds" % (end_time - start_time))
os.system('say "Your experiment has finished. Please collect your results"')
print(cnt)
