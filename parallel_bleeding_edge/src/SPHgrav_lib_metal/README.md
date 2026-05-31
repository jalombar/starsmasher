# SPHgrav_lib_metal

Direct-summation gravity library for Apple Silicon using Metal.

It exposes the same ABI used by the NVIDIA backend:

- `gpu_init_dev_`
- `firsthalf_grav_forces_`
- `lasthalf_grav_forces_`

The Apple Silicon GPU target is built from `parallel_bleeding_edge/src`:

```bash
make gpu
```

The resulting executable is written to:

- `parallel_bleeding_edge/bin/starsmasher_gpu`

The main Makefile selects this backend automatically on `Darwin/arm64`, so users should not need to edit any source files or Makefiles to build on Apple Silicon.

Do not force `DYLD_SHARED_REGION=private` when running the Metal backend. On macOS, that setting can prevent Metal from seeing the GPU device.
