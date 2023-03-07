# torchvision

PyTorch/1.11.0-foss-2021b-CUDA-11.4.1 is built with Pillow 8.3.2
```
>>> import PIL
>>> PIL.__version__
'8.3.2'
```
But torchvision requires anything but 8.3.x

```
 ('Pillow-SIMD', '8.2.0'), #  torchvision 0.13.0 has requirement pillow!=8.3.*,>=5.3.0, but you have pillow 8.3.2.
```
