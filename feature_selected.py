'''
Use Deepmg
Created: 21/02/2020
'''

import os
import sys
import time

if sys.platform == 'darwin':
    ROOT = '/Users/dcongtinh'
    path_data = ROOT + '/Downloads/data/'
else:
    ROOT = '/home/ubuntu'
    path_data = ROOT + '/data/'

dataset = ['cirgene', 'colgene', 'ibdgene', 'obegene', 't2dgene', 'wt2dgene']

# Set options as input for functions to reduce dim
algo_redu = 'pc'
new_dim = 576   # image 32x32 pixels
run_time = 10
num_bin = 10
type_bin = 'eqf'
scale_mode = 'qtf'
model_v = 'fc_model'

cnt = 0
def executeCommand(__command__):
    global cnt
    cnt += 1
    start_time = time.time()
    os.system(__command__)
    end_time = time.time()
    print("%f"%(end_time-start_time))
    # os.system('say Done %d' % cnt)
    # print(__command__ + '\n')

start_time = time.time()
time_start = time.ctime()

for dataset_name in dataset:
    # Set the PATH:
    file_run = ROOT +'/gene-abundance/subdeepmg/__main__.py' # Use package deepmg to modify
    file_read_results = ROOT +'/gene-abundance/subdeepmg/read_results.py' # Use package deepmg to modify
    machine_path = ROOT +'/gene-abundance' #for my macbook
    type_emb = 'fills'
    colormap = 'gray'

    res_folder = '%s/experiment/__results__/%s/%s_%s%d_10_%s%s_nb%d_auy_%s/' % (machine_path, model_v, scale_mode, algo_redu, new_dim, type_emb, type_bin, num_bin, colormap)
    img_folder = '%s/experiment/__images__/%s/%s_%s%d_10_%s%s_nb%d_auy_%s/' % (machine_path, model_v, scale_mode, algo_redu, new_dim, type_emb, type_bin, num_bin, colormap)

    command = 'python3 %s --save_para y --original_data_folder %s --data_name %s --run_time %d --algo_redu %s --new_dim %d --model %s --auto_v y --parent_folder_results %s --parent_folder_img %s --scale_mode %s --type_bin %s --num_bin %d --type_emb %s -z 255 --channel 1 --colormap %s --log_file %s' % (file_run, path_data, dataset_name, run_time, algo_redu, new_dim, model_v, res_folder, img_folder, scale_mode, type_bin, num_bin, type_emb, colormap, res_folder + 'feature_selected_' + dataset_name + '.txt')
    
    executeCommand(command)
    os.system('python3 %s -i %s -o %s' % (file_read_results, res_folder, res_folder))

print("Training dataset use [%s] at %s\n.\n.\n." % (algo_redu.upper(), time_start))
print("Run-time: %f seconds" % (end_time - start_time))
# os.system('say "Your experiment has finished. Please collect your results"')
print(cnt)
