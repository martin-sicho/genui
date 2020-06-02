"""
gpu

Created by: Martin Sicho
On: 6/2/20, 9:11 AM
"""

import nvgpu

ALLOCATIONS = set()

class GPUAllocationException(Exception):
    
    def __init__(device, *args, **kwargs):
        msg = f'Failed to allocate GPU device: {device}'
        super().__init__(msg, *args, **kwargs)

def check_availability():
    available = False
    try:
        nvgpu.gpu_info()
        available = True
    except FileNotFoundError as exp:
        if exp.filename == 'nvidia-smi':
            available = False
        else:
            raise exp

    return available

def info(device=None, memsort=False):
    if not check_availability():
        return []
    
    try:
        if type(device) == str:
            return next(x for x in nvgpu.gpu_info() if x['uuid'] == device)
        if type(device) == int:
            return next(x for x in nvgpu.gpu_info() if int(x['index']) == device)
    except StopIteration:
        raise Exception(f'Failed to find device: {device}')
    
    devices = nvgpu.gpu_info()
    if memsort:
        devices.sort( 
            key=lambda x: (1 - x['mem_used_percent'], x['mem_total']),
            reverse=True
        )
    return devices

def allocate():
    if not check_availability():
        return None

    device = info(memsort=True)[0]
    uuid = device['uuid']
    if uuid not in ALLOCATIONS:
        ALLOCATIONS.add(uuid)
    else:
        print(f'Device {uuid} is already allocated.')
        raise GPUAllocationException(uuid)
        
    return device

def release(device):
    try:
        ALLOCATIONS.remove(device['uuid'])
    except KeyError:
        print(f'Device already released: {device}')
