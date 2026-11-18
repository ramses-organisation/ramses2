## GPU performance ##

Measured speeds of the CUDA hydro and MHD solvers on one GPU with one MPI rank, to help choose the solver, the precision and `NSUBGRID` for a run. All numbers come from `develop` at 1fd9b705 (18 November 2026). Other problems, namelists and code versions will give different speeds.

Speeds are in millions of cell updates per second (Mcell/s) for the whole time step, not only the Godunov kernel. On AMR runs every updated cell is counted, on every level.

- "double" is `NPRE=8` with the default math.
- "single" is `NPRE=4` with `FASTMATH=1` (`-gpu=fastmath`).

The GPUs are an A100-PCIE-40GB (`CUDA_ARCH=sm_80`), an H200 NVL (`sm_90`) and a B200 (`sm_100`).

### Hydro (`MHD=0`, `riemann='hllc'`) ###

Uniform periodic 3D Sedov blast at 512^3, `slope_type=2`. Each value is the median coarse-step time over the last 16 of 32 steps, averaged over three runs.

| precision | `NSUBGRID` | A100 | H200 | B200 |
| --------- | ---------- | ---- | ---- | ---- |
| double | 1 | 739 | 1298 | 1746 |
| double | 2 | 1510 | 2757 | 3578 |
| double | 4 | does not fit | 3435 | 4456 |
| single | 1 | 1648 | 2159 | 3725 |
| single | 2 | 3374 | 6657 | 8040 |
| single | 4 | 3341 | 7994 | 10093 |

In double precision at `NSUBGRID=4` the kernel needs 207 KB of shared memory per block, more than the A100 has.

### MHD (`MHD=1`, `INIT=OT`) ###

3D Orszag-Tang vortex, `slope_type=2`, `slope_mag_type=2`, with `riemann2d='hlld'`. The `hlld` columns use `riemann='hlld'` and the `uct-hlld` columns use `riemann='uct-hlld'`. MHD runs with `NSUBGRID` 1 or 2.

**Uniform grid.** The grid is 512^3 on the H200 and the B200, and 256^3 on the A100, where `uct-hlld` at 512^3 does not fit in double precision. Each value is the median coarse-step time over the last 16 of 32 steps, averaged over three runs.

| precision | `NSUBGRID` | A100 `hlld` | A100 `uct-hlld` | H200 `hlld` | H200 `uct-hlld` | B200 `hlld` | B200 `uct-hlld` |
| --------- | ---------- | ----------- | --------------- | ----------- | --------------- | ----------- | --------------- |
| double | 1 | 138 | 237 | 266 | 393 | 310 | 508 |
| double | 2 | 241 | 465 | 400 | 869 | 530 | 860 |
| single | 1 | 378 | 518 | 667 | 804 | 828 | 1155 |
| single | 2 | 791 | 984 | 1482 | 2096 | 1772 | 2473 |

**AMR.** Levels 7 to 9, refined on the field gradient (`err_grad_B=0.15`), with `nsubcycle=1` on every level. Levels 8 and 9 cover about 35% and 9% of the volume. Each value is the median coarse-step time over steps 41 to 80, averaged over two runs. The runs use `ngridmax=4000000` and `ncachemax=1000000`, and `ngridmax=3400000` and `ncachemax=600000` on the A100.

| precision | `NSUBGRID` | A100 `hlld` | A100 `uct-hlld` | H200 `hlld` | H200 `uct-hlld` | B200 `hlld` | B200 `uct-hlld` |
| --------- | ---------- | ----------- | --------------- | ----------- | --------------- | ----------- | --------------- |
| double | 1 | 96 | 138 | 168 | 227 | 198 | 296 |
| double | 2 | 156 | 214 | 215 | 409 | 268 | 461 |
| single | 1 | 208 | 235 | 385 | 398 | 500 | 548 |
| single | 2 | 348 | 369 | 610 | 669 | 786 | 840 |

### GPU memory ###

Peak GPU memory of the 512^3 uniform runs on the H200, as sampled by `nvidia-smi`:

| run | double | single |
| --- | ------ | ------ |
| hydro, `NSUBGRID=1` | 15.1 GiB | 9.4 GiB |
| MHD `hlld`, `NSUBGRID=1` | 28.9 GiB | 16.3 GiB |
| MHD `uct-hlld`, `NSUBGRID=1` | 66.6 GiB | 35.2 GiB |

`uct-hlld` is faster because it solves each face once and keeps the results in a GPU buffer, which `hlld` does not need. The buffer holds about 2 KB per oct slot in double precision (1 KB in single) for all `ngridmax + ncachemax` slots, whether or not they are in use. If it cannot be allocated, the run stops and prints the size it needs.
