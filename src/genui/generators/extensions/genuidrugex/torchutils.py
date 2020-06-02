"""
torchutils

Created by: Martin Sicho
On: 6/2/20, 9:13 AM
"""

import torch

def pretty_size(size):
	"""
	Pretty prints a torch.Size object

	By user machinethink: https://forums.fast.ai/t/gpu-memory-not-being-freed-after-training-is-over/10265/7

	"""
	assert(isinstance(size, torch.Size))
	return " × ".join(map(str, size))

def dump_tensors(gpu_only=True):
	"""
	Prints a list of the Tensors being tracked by the garbage collector.

	By user machinethink: https://forums.fast.ai/t/gpu-memory-not-being-freed-after-training-is-over/10265/7

	"""
	import gc
	total_size = 0
	for obj in gc.get_objects():
		try:
			if torch.is_tensor(obj):
				if not gpu_only or obj.is_cuda:
					print("%s:%s%s %s" % (type(obj).__name__,
										  " GPU" if obj.is_cuda else "",
										  " pinned" if obj.is_pinned else "",
										  pretty_size(obj.size())))
					total_size += obj.numel()
			elif hasattr(obj, "data") and torch.is_tensor(obj.data):
				if not gpu_only or obj.is_cuda:
					print("%s → %s:%s%s%s%s %s" % (type(obj).__name__,
												   type(obj.data).__name__,
												   " GPU" if obj.is_cuda else "",
												   " pinned" if obj.data.is_pinned else "",
												   " grad" if obj.requires_grad else "",
												   " volatile" if obj.volatile else "",
												   pretty_size(obj.data.size())))
					total_size += obj.data.numel()
		except Exception as e:
			pass
	print("Total size cleaned:", total_size)

def cleanup():
	dump_tensors()
	torch.cuda.empty_cache()