"""
gpu

Created by: Martin Sicho
On: 6/2/20, 9:11 AM
"""

import nvgpu

def check_gpu_availability():
    available = False
    try:
        nvgpu.gpu_info()
    except FileNotFoundError as exp:
        if exp.filename == 'nvidia-smi':
            available = False
        else:
            raise exp

    return available


