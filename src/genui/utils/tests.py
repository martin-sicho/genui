from django.test import TestCase
from unittest import skipIf
from . import gpu

class GPUTests(TestCase):

    @skipIf(not gpu.check_availability(), "No CUDA enabled gpus were found on the system. Skipping test...")
    def test_info(self):
        device_info = gpu.info(device=0)
        self.assertTrue(device_info['index'] == '0')
        devices = gpu.info(memsort=True)
        print(devices)
        self.assertTrue(max([x['mem_used_percent'] for x in devices]) == devices[0]['mem_used_percent'])
