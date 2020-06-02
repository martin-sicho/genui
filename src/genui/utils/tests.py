from django.test import TestCase
from unittest import skipIf
from . import gpu

class GPUTests(TestCase):

    @skipIf(not gpu.check_gpu_availability(), "No CUDA enabled gpus were found on the system. Skipping test...")
    def test_info(self):
        # TODO: implement
        pass
