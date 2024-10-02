import numpy as np
import subprocess
import os
import time

def write_arr(arr, f):
    np.savetxt(f, arr, fmt='%f')

def write_arr_int(arr, f):
    np.savetxt(f, arr, fmt='%d')

def write_kernel(f, cpp_ker, sigma_w):
    f.write(f"Type={cpp_ker['Type']}\n")
    f.write(f"Function={cpp_ker['Function']}\n")
    f.write(f"SigmaW={sigma_w}\n")

def read_result(file_results):
    # Implement the logic to read results from the file
    pass

def match_cpp(targets, s, file_results, file_temp, match_float_cmd, match_double_cmd, num_rand):
    with open(file_temp, 'w') as f:
        for target in targets:
            if target['type'] == 'curvacc':
                f.write('CurveAcc\n')
                f.write(f"Range=\n{np.min(target['vx'])} {np.max(target['vx'])}\n")
                f.write(f"Weight=\n{target['weight']}\n")
                if 'CppKer' not in target:
                    target['CppKer'] = {'Type': 'SqDistScalar', 'Function': 'Gaussian'}
                    print(f"no CppKer field given for target {i}. Assuming scalar gaussian kernel")
                f.write('Kernel=\n')
                write_kernel(f, target['CppKer'], target['sigmaW'])
                if target['CppKer']['Type'] in ['CauchyGpu', 'GaussGpu'] and s.get('typefloat') != 'float':
                    s['typefloat'] = 'float'
                    print("redefining typefloat to 'float' for use of Gpu code")
                f.write('VY=\n')
                write_arr(target['y'], f)
                f.write('FX=\n')
                write_arr_int(target['vx'] - np.min(target['vx']) + 1, f)
                f.write('FY=\n')
                write_arr_int(target['vy'], f)
                f.write('WX,WY=\n')
                if 'wx' not in target:
                    target['wx'] = np.ones(target['vx'].shape[1])
                write_arr(target['wx'], f)
                if 'wy' not in target:
                    target['wy'] = np.ones(target['vy'].shape[1])
                write_arr(target['wy'], f)
            elif target['type'] == 'landmarks':
                f.write('Landmarks\n')
                f.write(f"Range=\n{target['vx'][0]} {target['vx'][-1]}\n")
                f.write(f"Weight=\n{target['weight']}\n")
                f.write('Y=\n')
                write_arr(target['y'], f)
            elif target['type'] == 'l2image':
                f.write('L2Image\n')
                f.write(f"Range=\n{target['rx'][0]} {target['rx'][-1]}\n")
                f.write(f"Weight=\n{target['weight']}\n")
                f.write('SourceImage=\n')
                write_arr(target['imsource'], f)
                f.write('TargetImage=\n')
                write_arr(target['imtarget'], f)
                f.write(f"TargetGridBase=\n{target['basetarget'][0]} {target['basetarget'][1]} {target['basetarget'][2]}\n")
                f.write(f"TargetGridVoxSize=\n{target['voxsizetarget'][0]} {target['voxsizetarget'][1]} {target['voxsizetarget'][2]}\n")

    s.setdefault('optim_maxiter', 500)
    s.setdefault('optim_stepsize', 1)
    s.setdefault('optim_breakratio', 1e-6)
    s.setdefault('optim_loopbreak', 40)
    s.setdefault('typefloat', 'double')

    if os.path.exists(file_results):
        os.remove(file_results)

    match_cmd = match_float_cmd if s['typefloat'] == 'float' else match_double_cmd

    s.setdefault('optim_useoptim', 'adaptdesc')
    useoptim = 1 if s['optim_useoptim'] == 'adaptdesc' else 0

    command = [
        match_cmd, '-d', file_temp, '-o', file_results, '-i', str(s['optim_maxiter']),
        '-s', str(s['optim_stepsize']), '-w', str(s['gammaR']), '-u', str(useoptim),
        '-b', str(s['optim_breakratio']), '-l', str(s['optim_loopbreak'])
    ]

    start_time = time.time()
    subprocess.run(command)
    elapsed_time = time.time() - start_time

    s['typefloat'] = s['typefloat']
    s.update(read_result(file_results))
    s['elapsedTime'] = elapsed_time

    for file in os.listdir('.'):
        if file.startswith(num_rand):
            os.remove(file)


def blobimage():
    return None